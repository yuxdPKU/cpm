./run_cpm_b_chain.sh \
  --input run79516list.txt \
  --input-is-list \
  --out-dir output/run79516 \
  --prefix run79516 \
  --b2-containers \
  --b2-unweighted \
  --no-combined-output \
  --no-keep-intermediates

./run_cpm_b_chain.sh \
  --input run79516list.txt \
  --input-is-list \
  --out-dir output/run79516 \
  --prefix run79516_ptweight \
  --b2-containers \
  --b2-weighted \
  --no-combined-output \
  --no-keep-intermediates

./run_cpm_b_chain.sh \
  --input run79516list.txt \
  --input-is-list \
  --out-dir output/run79516 \
  --prefix run79516_qaoffline \
  --b1-max-pair-records 0 \
  --no-b1-write-pairs \
  --b2-input-mode batches \
  --b2-unweighted \
  --no-combined-output \
  --keep-intermediates

./run_cpm_b_chain.sh \
  --input run79516list.txt \
  --input-is-list \
  --out-dir output/run79516 \
  --prefix run79516_ptweight_qaoffline \
  --b1-max-pair-records 0 \
  --no-b1-write-pairs \
  --b2-input-mode batches \
  --b2-weighted \
  --no-combined-output \
  --keep-intermediates
