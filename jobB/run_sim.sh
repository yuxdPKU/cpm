./run_cpm_b_chain.sh \
  --input simlist.txt \
  --input-is-list \
  --out-dir output/sim \
  --prefix sim \
  --b2-containers \
  --b2-unweighted \
  --no-combined-output \
  --no-keep-intermediates

./run_cpm_b_chain.sh \
  --input simlist.txt \
  --input-is-list \
  --out-dir output/sim \
  --prefix sim_ptweight \
  --b2-containers \
  --b2-weighted \
  --no-combined-output \
  --no-keep-intermediates

./run_cpm_b_chain.sh \
  --input simlist.txt \
  --input-is-list \
  --out-dir output/sim \
  --prefix sim_qaoffline \
  --b1-max-pair-records 0 \
  --no-b1-write-pairs \
  --b2-input-mode batches \
  --b2-unweighted \
  --no-combined-output \
  --keep-intermediates

./run_cpm_b_chain.sh \
  --input simlist.txt \
  --input-is-list \
  --out-dir output/sim \
  --prefix sim_qaoffline_ptweight \
  --b1-max-pair-records 0 \
  --no-b1-write-pairs \
  --b2-input-mode batches \
  --b2-weighted \
  --no-combined-output \
  --keep-intermediates
