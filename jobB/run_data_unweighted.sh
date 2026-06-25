#!/usr/bin/env bash
set -euo pipefail

echo run data unweighted
./run_cpm_b_chain.sh \
  --input run79516list.txt \
  --input-is-list \
  --out-dir output/run79516_unweighted \
  --prefix run79516_unweighted_maxdca0p2 \
  --max-pair-dca 0.2 \
  --unweighted

#echo run data qa unweighted
#./run_cpm_qa_chain.sh \
#  --input run79516list.txt \
#  --input-is-list \
#  --out-dir output/run79516_unweighted_qa \
#  --prefix run79516_unweighted \
#  --unweighted \
#  --write-pair-tree
