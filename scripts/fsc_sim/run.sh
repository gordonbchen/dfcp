#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/../.." && pwd)
data_dir="${1:-data/fsc}"
output_dir="${2:-$data_dir/prepared}"
fsc="${FSC:-fsc28}"
threads="${THREADS:-8}"
seed="${SEED:-756789}"
template="${TEMPLATE:-ex_0_pop_1.tpl}"
est="${EST:-ex_0_pop.est}"
simulation=${template%.tpl}

cd "$repo_root"
(
    cd "$data_dir"
    "$fsc" -t "$template" -n1 -e "$est" -T -E1 -s0 -c"$threads" -I -G -x --seed "$seed"
)

python3 "$script_dir/prep_data.py" \
    "$data_dir/$simulation/${simulation}_1_1.gen" \
    "$output_dir"
