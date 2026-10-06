runs=( 79516 )
types=(
	'CPM_data_reco_FITacts_USEMMStrue'
	'CPM_data_reco_FITacts_USEMMSfalse'
)
#runs=( 29 )
#types=(
#	'CPM_sim_reco_FITacts_DISTORTIONINPUTtrue_TRUTHSEEDINGtrue'
#	'CPM_sim_reco_FITacts_DISTORTIONINPUTtrue_TRUTHSEEDINGfalse'
#	'CPM_sim_reco_FITacts_DISTORTIONINPUTfalse_TRUTHSEEDINGtrue'
#	'CPM_sim_reco_FITacts_DISTORTIONINPUTfalse_TRUTHSEEDINGfalse'
#	'CPM_sim_reco_FITgenfit_DISTORTIONINPUTtrue_TRUTHSEEDINGtrue'
#	'CPM_sim_reco_FITgenfit_DISTORTIONINPUTtrue_TRUTHSEEDINGfalse'
#	'CPM_sim_reco_FITgenfit_DISTORTIONINPUTfalse_TRUTHSEEDINGtrue'
#	'CPM_sim_reco_FITgenfit_DISTORTIONINPUTfalse_TRUTHSEEDINGfalse'
#	'CPM_sim_reco_FITtruth_DISTORTIONINPUTtrue_TRUTHSEEDINGtrue'
#	'CPM_sim_reco_FITtruth_DISTORTIONINPUTtrue_TRUTHSEEDINGfalse'
#	'CPM_sim_reco_FITtruth_DISTORTIONINPUTfalse_TRUTHSEEDINGtrue'
#	'CPM_sim_reco_FITtruth_DISTORTIONINPUTfalse_TRUTHSEEDINGfalse'
#)
echo nruns ${#runs[@]}
echo ntypes ${#types[@]}
for ((k=0; k<${#runs[@]}; k++))
do
  for ((j=0; j<${#types[@]}; j++))
  do
    echo run ${runs[$k]} ${types[$j]}
    find /sphenix/u/xyu3/workarea/cpm/root/Reconstructed/${runs[$k]} \
      -maxdepth 1 \
      -name "${types[$j]}_${runs[$k]}-[0-9]*.root_CPMVoxelContainer.root" \
      -print | sort -V > ./list/list_run${runs[$k]}_${types[$j]}.txt
  done
done
#
#find /sphenix/u/xyu3/workarea/cpm/root/Reconstructed/29 \
#  -maxdepth 1 \
#  -name 'CPM_sim_reco_acts_29-*.root_CPMVoxelContainer.root' \
#  -print | sort -V > list_CPM_run29_sim_acts.txt
#
#find /sphenix/u/xyu3/workarea/cpm/root/Reconstructed/29 \
#  -maxdepth 1 \
#  -name 'CPM_sim_reco_genfit_29-*.root_CPMVoxelContainer.root' \
#  -print | sort -V > list_CPM_run29_sim_genfit.txt
#
#find /sphenix/u/xyu3/workarea/cpm/root/Reconstructed/29 \
#  -maxdepth 1 \
#  -name 'CPM_sim_reco_truth_29-*.root_CPMVoxelContainer.root' \
#  -print | sort -V > list_CPM_run29_sim_truth.txt
#
#find /sphenix/u/xyu3/workarea/cpm/root/Reconstructed/29 \
#  -maxdepth 1 \
#  -name 'CPM_sim_reco_truth_withdistortion_29-*.root_CPMVoxelContainer.root' \
#  -print | sort -V > list_CPM_run29_sim_truth_withdistortion.txt
#
#
#if list has already generated, use following command
#sort -V -o run79516list.txt run79516list.txt
