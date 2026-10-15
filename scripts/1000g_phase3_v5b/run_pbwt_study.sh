#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
windows_dir=${1:-data/1000g_phase3_v5b/windows_8_maf1}
output_dir=${2:-output/maf1_pbwt_study}
n_parallel=${N_PARALLEL:-4}

cd "$repo_root"
mkdir -p "$output_dir"

for k in 128 256 512 1024 2048; do
    name="maf1_joint_k${k}_step1"
    save_assign=()
    if (( k == 1024 )); then
        save_assign=(--save-r-assign)
    fi
    "$script_dir/impute.sh" "$windows_dir" \
        --name "$name" --n-parallel "$n_parallel" "${save_assign[@]}" \
        --init emission --target_fit joint --block_max_k "$k" --max_train_steps 1
    ./build/eval_impute "$windows_dir" "$name" "$output_dir/$name.tsv"
done

for steps in 0 3 30; do
    name="maf1_joint_k1024_step${steps}"
    if (( steps == 30 )); then
        name=maf1_joint_k1024_converged
    fi
    "$script_dir/impute.sh" "$windows_dir" \
        --name "$name" --n-parallel "$n_parallel" --save-r-assign \
        --init emission --target_fit joint --block_max_k 1024 --max_train_steps "$steps"
    ./build/eval_impute "$windows_dir" "$name" "$output_dir/$name.tsv"
done
