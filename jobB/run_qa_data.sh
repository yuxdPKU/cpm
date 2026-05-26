./run_cpm_qa_chain.sh \
  --input run79516list.txt \
  --input-is-list \
  --out-dir output/run79516_qa \
  --prefix run79516 \
  --unweighted

./run_cpm_qa_chain.sh \
  --input run79516list.txt \
  --input-is-list \
  --out-dir output/run79516_qa \
  --prefix run79516_ptweight \
  --weighted
