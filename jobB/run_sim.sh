./run_cpm_b_chain.sh \
  --input simlist.txt \
  --input-is-list \
  --out-dir output/sim \
  --prefix sim \
  --unweighted

./run_cpm_b_chain.sh \
  --input simlist.txt \
  --input-is-list \
  --out-dir output/sim \
  --prefix sim_ptweight \
  --weighted

./run_cpm_b_chain.sh \
  --input simlist.txt \
  --input-is-list \
  --out-dir output/sim \
  --prefix sim_qaoffline \
  --unweighted

./run_cpm_b_chain.sh \
  --input simlist.txt \
  --input-is-list \
  --out-dir output/sim \
  --prefix sim_qaoffline_ptweight \
  --weighted
