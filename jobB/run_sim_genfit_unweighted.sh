#!/usr/bin/env bash
set -euo pipefail

echo run sim genfit unweighted
./run_cpm_b_chain.sh \
  --input list_CPM_run29_sim_genfit.txt \
  --input-is-list \
  --out-dir output/sim_run29_genfit_unweighted \
  --prefix sim_run29_genfit_unweighted_maxdca0p2 \
  --max-pair-dca 0.2 \
  --unweighted

#echo run sim genfit qa unweighted
#./run_cpm_qa_chain.sh \
#  --input simlist_genfit.txt \
#  --input-is-list \
#  --out-dir output/sim_genfit_unweighted_qa \
#  --prefix sim_genfit_unweighted \
#  --unweighted \
#  --write-pair-tree
