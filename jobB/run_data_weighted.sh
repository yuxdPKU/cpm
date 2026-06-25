#!/usr/bin/env bash
set -euo pipefail

echo run data weighted
./run_cpm_b_chain.sh \
  --input run79516list.txt \
  --input-is-list \
  --out-dir output/run79516_weighted \
  --prefix run79516_weighted_maxdca0p2 \
  --max-pair-dca 0.2 \
  --weighted

#echo run data qa weighted
#./run_cpm_qa_chain.sh \
#  --input run79516list.txt \
#  --input-is-list \
#  --out-dir output/run79516_weighted_qa \
#  --prefix run79516_weighted \
#  --weighted \
#  --write-pair-tree
