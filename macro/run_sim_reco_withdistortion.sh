#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n new  # setup sPHENIX environment in the singularity container shell. Note the shell is bash by default

# Additional commands for my local environment
export SPHENIX=/sphenix/u/xyu3
export MYINSTALL=$SPHENIX/install

# Setup MYINSTALL to local directory and run sPHENIX setup local script
# to adjust PATH, LD LIBRARY PATH, ROOT INCLUDE PATH, etc
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

echo "sPHENIX environment setup finished"

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`
echo rsyncing from $this_dir
echo running: $this_script $*

if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]
then
  cd $_CONDOR_SCRATCH_DIR
  rsync -av $this_dir/* .
else
  echo condor scratch NOT set
  exit -1
fi

nEvents=$1
InDst=$2
OutDir=$3
OutPrefix=$4
Index=$5
StepSize=$6
FitMode=${7:-auto}
UseTruthInformation=${8:-false}

#getinputfiles.pl $InDst $InSeedDst
getinputfiles.pl --filelist $InDstList

# print the environment - needed for debugging
printenv

root.exe -q -b Fun4All_G4_sPHENIX_reco_withdistortion.C\($nEvents,\"${InDst}\",\"${OutDir}\",\"${OutPrefix}\",$Index,$StepSize,true,false,\"${FitMode}\",${UseTruthInformation}\)
echo Script done
