#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
windows_dir=${1:?Usage: impute.sh WINDOWS_DIR [OPTION VALUE]...}
shift

cd "$repo_root"
./build.sh

shopt -s nullglob
window_dirs=("$windows_dir"/window_*)
if (( ${#window_dirs[@]} == 0 )); then
    echo "error: no window directories in $windows_dir" >&2
    exit 1
fi

for window_dir in "${window_dirs[@]}"; do
    echo "imputing: $window_dir" >&2
    ./build/impute \
        "$window_dir/ref.bin" \
        "$window_dir/target_observed.bin" \
        "$window_dir/observed_loci.txt" \
        "$window_dir/probs.bin" \
        "$@" > "$window_dir/impute.json"
done
