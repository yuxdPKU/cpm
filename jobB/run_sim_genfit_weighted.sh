#!/usr/bin/env bash
set -euo pipefail

echo run sim genfit weighted
./run_cpm_b_chain.sh \
  --input simlist_genfit.txt \
  --input-is-list \
  --out-dir output/sim_weighted \
  --prefix sim_genfit_weighted \
  --weighted

echo run sim genfit qa weighted
./run_cpm_qa_chain.sh \
  --input simlist_genfit.txt \
  --input-is-list \
  --out-dir output/sim_weighted_qa \
  --prefix sim_genfit_weighted \
  --weighted \
  --write-pair-tree
