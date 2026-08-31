#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
windows_dir=${1:?Usage: impute.sh WINDOWS_DIR [OPTION VALUE]...}
shift
n_parallel=1
impute_args=()
while (( $# )); do
    if [[ $1 == --n-parallel ]]; then
        n_parallel=${2:?--n-parallel requires a value}
        shift 2
    else
        impute_args+=("$1")
        shift
    fi
done

cd "$repo_root"
./build.sh

shopt -s nullglob
window_dirs=("$windows_dir"/window_*)
if (( ${#window_dirs[@]} == 0 )); then
    echo "error: no window directories in $windows_dir" >&2
    exit 1
fi

run_window() {
    local window_dir=$1
    shift
    exec ./build/impute \
        "$window_dir/ref.bin" \
        "$window_dir/target_observed.bin" \
        "$window_dir/observed_loci.txt" \
        "$window_dir/probs.bin" \
        "$@" > "$window_dir/impute.json"
}

stop_jobs() {
    trap - INT TERM
    local pids
    mapfile -t pids < <(jobs -pr)
    if (( ${#pids[@]} )); then
        kill "${pids[@]}" 2>/dev/null || true
        wait "${pids[@]}" 2>/dev/null || true
    fi
    exit 130
}
trap stop_jobs INT TERM

active=0
failed=0
declare -A job_dirs
declare -A job_logs
wait_one() {
    local pid
    if ! wait -n -p pid; then failed=1; fi
    printf '\n=== %s ===\n' "${job_dirs[$pid]}" >&2
    cat "${job_logs[$pid]}" >&2
    unset 'job_dirs[$pid]' 'job_logs[$pid]'
    active=$((active - 1))
}

for window_dir in "${window_dirs[@]}"; do
    log_file="$window_dir/impute.log"
    run_window "$window_dir" "${impute_args[@]}" 2> "$log_file" &
    pid=$!
    job_dirs[$pid]=$window_dir
    job_logs[$pid]=$log_file
    active=$((active + 1))
    if (( active == n_parallel )); then wait_one; fi
done
while (( active )); do wait_one; done
exit "$failed"
