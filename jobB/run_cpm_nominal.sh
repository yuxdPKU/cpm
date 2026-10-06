#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n new  # setup sPHENIX environment in the singularity container shell. Note the shell is bash by default

# Additional commands for my local environment
export SPHENIX=/sphenix/u/xyu3
export MYINSTALL=$SPHENIX/install

# Setup MYINSTALL to local directory and run sPHENIX setup local script
# to adjust PATH, LD LIBRARY PATH, ROOT INCLUDE PATH, etc
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

echo "sPHENIX environment setup finished"

INPUT=$1
IS_LIST=$2
OUT_DIR=$3
OUT_PREFIX=$4
DO_WEIGHT=$5
MAX_DCA=$6
DO_QA=$7

WEIGHT_OPTION="--unweighted"
if [[ "$DO_WEIGHT" == "true" ]]; then
  WEIGHT_OPTION="--weighted"
fi

LIST_OPTION=""
if [[ "$IS_LIST" == "true" ]]; then
  LIST_OPTION="--input-is-list"
fi

./run_cpm_b_chain.sh \
  --input "$INPUT" \
  $LIST_OPTION \
  --out-dir "$OUT_DIR" \
  --prefix "$OUT_PREFIX" \
  --max-pair-dca "$MAX_DCA" \
  $WEIGHT_OPTION

if [[ "$DO_QA" == "true" ]]; then
  ./run_cpm_qa_chain.sh \
    --input "$INPUT" \
    $LIST_OPTION \
    --out-dir "${OUT_DIR}_qa" \
    --prefix "$OUT_PREFIX" \
    --max-pair-dca "$MAX_DCA" \
    $WEIGHT_OPTION \
    --write-pair-tree
fi
