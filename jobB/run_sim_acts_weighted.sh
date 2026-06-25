#!/usr/bin/env bash
set -euo pipefail

echo run sim acts weighted
./run_cpm_b_chain.sh \
  --input list_CPM_run29_sim_acts.txt \
  --input-is-list \
  --out-dir output/sim_run29_acts_weighted \
  --prefix sim_run29_acts_weighted_maxdca0p2 \
  --max-pair-dca 0.2 \
  --weighted

#echo run sim acts qa weighted
#./run_cpm_qa_chain.sh \
#  --input simlist_acts.txt \
#  --input-is-list \
#  --out-dir output/sim_acts_weighted_qa \
#  --prefix sim_acts_weighted \
#  --weighted \
#  --write-pair-tree
