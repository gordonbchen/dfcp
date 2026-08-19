#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
threads="${THREADS:-1}"
n_windows="${N_WINDOWS:-200}"
overlap="${OVERLAP:-32}"
n_generate="${N_GENERATE:-$n_windows}"
windows_dir="${WINDOWS_DIR:-data/1000g_phase3_v5b/windows}"

cd "$repo_root"

"$script_dir/download.sh" "$@"
"$script_dir/biallelic_snvs.sh"
python3 "$script_dir/split_ref_target.py" --threads "$threads"
"$script_dir/mask_target.sh"
python3 "$script_dir/window.py" \
    --output-dir "$windows_dir" \
    --n-windows "$n_windows" \
    --overlap "$overlap" \
    --n-generate "$n_generate" \
    --threads "$threads"
"$script_dir/dfcp_prep.sh" "$windows_dir"
