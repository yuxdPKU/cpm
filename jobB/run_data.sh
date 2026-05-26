./run_cpm_b_chain.sh \
  --input run79516list.txt \
  --input-is-list \
  --out-dir output/run79516 \
  --prefix run79516 \
  --unweighted

./run_cpm_b_chain.sh \
  --input run79516list.txt \
  --input-is-list \
  --out-dir output/run79516 \
  --prefix run79516_ptweight \
  --weighted

./run_cpm_b_chain.sh \
  --input run79516list.txt \
  --input-is-list \
  --out-dir output/run79516 \
  --prefix run79516_qaoffline \
  --unweighted

./run_cpm_b_chain.sh \
  --input run79516list.txt \
  --input-is-list \
  --out-dir output/run79516 \
  --prefix run79516_ptweight_qaoffline \
  --weighted
