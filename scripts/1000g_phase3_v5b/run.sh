#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
threads="${THREADS:-1}"

cd "$repo_root"

"$script_dir/download.sh" "$@"
"$script_dir/biallelic_snvs.sh"
python3 "$script_dir/split_ref_target.py" --threads "$threads"
"$script_dir/mask_target.sh"
"$script_dir/dfcp_prep.sh"
