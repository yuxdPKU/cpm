#!/usr/bin/env bash
set -euo pipefail

echo run sim acts unweighted
./run_cpm_b_chain.sh \
  --input list_CPM_run29_sim_acts.txt \
  --input-is-list \
  --out-dir output/sim_run29_acts_unweighted \
  --prefix sim_run29_acts_unweighted_maxdca0p2 \
  --max-pair-dca 0.2 \
  --unweighted

#echo run sim acts qa unweighted
#./run_cpm_qa_chain.sh \
#  --input simlist_acts.txt \
#  --input-is-list \
#  --out-dir output/sim_acts_unweighted_qa \
#  --prefix sim_acts_unweighted \
#  --unweighted \
#  --write-pair-tree
