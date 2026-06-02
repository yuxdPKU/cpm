#!/usr/bin/env bash
set -euo pipefail

echo run sim genfit unweighted
./run_cpm_b_chain.sh \
  --input simlist_genfit.txt \
  --input-is-list \
  --out-dir output/sim_genfit_unweighted \
  --prefix sim_genfit_unweighted \
  --unweighted

echo $?

echo run sim genfit qa unweighted
./run_cpm_qa_chain.sh \
  --input simlist_genfit.txt \
  --input-is-list \
  --out-dir output/sim_genfit_unweighted_qa \
  --prefix sim_genfit_unweighted \
  --unweighted

echo $?
