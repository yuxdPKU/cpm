#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n new  # setup sPHENIX environment in the singularity container shell. Note the shell is bash by default

# Additional commands for my local environment
export SPHENIX=/sphenix/u/xyu3
export MYINSTALL=$SPHENIX/install

# Setup MYINSTALL to local directory and run sPHENIX setup local script
# to adjust PATH, LD LIBRARY PATH, ROOT INCLUDE PATH, etc
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

echo "sPHENIX environment setup finished"

useScratch=true

if [[ "${useScratch}" == true ]]; then
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
fi

nEvents=$1
InDst=$2
OutDir=$3
OutPrefix=$4
FitMode=$5
UseMMS=$6
WriteMiniDst=$7
WritePrunedSeedsToMiniDst=$8
Index=$9
StepSize=${10}

if [[ "${useScratch}" == true ]]; then
  if [[ "${InDst}" == *.root ]]; then
    getinputfiles.pl $InDst
  elif [[ "${InDst}" == *.list ]]; then
    getinputfiles.pl --filelist $InDst
  fi
fi

#if ! getinputfiles.pl "$InDst"
#then
#  echo "ERROR: getinputfiles.pl failed" >&2
#  exit 1
#fi
#
#if [[ ! -s "$InDst" ]]
#then
#  echo "ERROR: missing seed DST: $InDst" >&2
#  exit 1
#fi

# print the environment - needed for debugging
printenv

if ! root.exe -q -b Fun4All_TrackAnalysis_CPM.C\($nEvents,\"${InDst}\",\"${OutDir}\",\"${OutPrefix}\",\"${FitMode}\",${UseMMS},${WriteMiniDst},${WritePrunedSeedsToMiniDst},${Index},${StepSize}\)
then
  echo "ERROR: ROOT reconstruction failed" >&2
  exit 1
fi
echo Script done
