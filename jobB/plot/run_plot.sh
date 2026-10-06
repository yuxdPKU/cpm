phiType=( 'central' 'west' 'east' )
types=( 'acts' 'genfit' 'truth' )
echo ${#phiType[@]}
for ((k=0; k<${#phiType[@]}; k++))
do
  echo run ${phiType[$k]}
  root -b -q -l draw1D_r_from3D_sim_closure_unweighted.C"(\"${phiType[$k]}\")"
  root -b -q -l draw1D_r_from3D_sim_closure_weighted.C"(\"${phiType[$k]}\")"
  root -b -q -l draw1D_z_from3D_sim_closure_unweighted.C"(\"${phiType[$k]}\")"
  root -b -q -l draw1D_z_from3D_sim_closure_weighted.C"(\"${phiType[$k]}\")"
  for ((j=0; j<${#types[@]}; j++))
  do
    echo run ${phiType[$k]} ${types[$j]}
    root -b -q -l draw2D_rz_from3D_sim_closure_unweighted.C"(29,\"${types[$j]}\",\"${phiType[$k]}\")"
    root -b -q -l draw2D_rz_from3D_sim_closure_weighted.C"(29,\"${types[$j]}\",\"${phiType[$k]}\")"
  done
done
