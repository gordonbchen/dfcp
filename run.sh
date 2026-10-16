#!/usr/bin/env bash
set -euo pipefail

./build.sh
printf '\n\n'

./build/impute data/fsc/prepared/ref.bin --output_r_assign output/r_assign.bin \
  --pbwt_match_len 20 \
  > /dev/null
printf '\n\n'

./build/eval_clusters data/fsc/prepared/ref.bin \
  output/r_assign.bin data/fsc/prepared/variant_pos.txt \
  data/fsc/ex_0_pop_1/ex_0_pop_1_1_true_trees.trees --tree_vis output/tree.dot \
  --tract_lens output/tract_lens.txt --clade_times output/clade_times.txt > /dev/null
printf '\n\n'

dot -Tsvg output/tree.dot -o output/tree20.svg
