./run_cpm_qa_chain.sh \
  --input simlist.txt \
  --input-is-list \
  --out-dir output/sim_qa \
  --prefix sim \
  --unweighted

./run_cpm_qa_chain.sh \
  --input simlist.txt \
  --input-is-list \
  --out-dir output/sim_qa \
  --prefix sim_ptweight \
  --weighted
