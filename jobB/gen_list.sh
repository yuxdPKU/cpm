find /sphenix/u/xyu3/workarea/cpm/root/Reconstructed/29 \
  -maxdepth 1 \
  -name 'CPM_sim_reco_acts_29-*.root_CPMVoxelContainer.root' \
  -print | sort -V > list_CPM_run29_sim_acts.txt

find /sphenix/u/xyu3/workarea/cpm/root/Reconstructed/29 \
  -maxdepth 1 \
  -name 'CPM_sim_reco_genfit_29-*.root_CPMVoxelContainer.root' \
  -print | sort -V > list_CPM_run29_sim_genfit.txt

find /sphenix/u/xyu3/workarea/cpm/root/Reconstructed/29 \
  -maxdepth 1 \
  -name 'CPM_sim_reco_truth_29-*.root_CPMVoxelContainer.root' \
  -print | sort -V > list_CPM_run29_sim_truth.txt

find /sphenix/u/xyu3/workarea/cpm/root/Reconstructed/29 \
  -maxdepth 1 \
  -name 'CPM_sim_reco_truth_withdistortion_29-*.root_CPMVoxelContainer.root' \
  -print | sort -V > list_CPM_run29_sim_truth_withdistortion.txt


#if list has already generated, use following command
#sort -V -o run79516list.txt run79516list.txt
