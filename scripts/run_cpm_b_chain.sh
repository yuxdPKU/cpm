#!/usr/bin/env bash

if [ -z "${BASH_VERSION:-}" ]; then
  exec /usr/bin/env bash "$0" "$@"
fi

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/run_cpm_b_chain.sh --input JOB_A.root [options]
  scripts/run_cpm_b_chain.sh --input cpm_filelist.txt --input-is-list [options]

Runs the CPM Job B macro chain:
  optional B0 build/check event index QA
  B1 local line-line PoCA
  B2 weighted voxel accumulator
  B3 average-correction histogram writer
  B3 histogram check
  combined Job B ROOT file

Options:
  --input PATH                  Job A CPMVoxelContainer ROOT file, or file list.
  --input-is-list               Treat --input as a newline-separated file list.
  --out-dir DIR                 Output directory. Default: .
  --prefix NAME                 Output filename prefix. Default: CPM
  --metadata PATH               Job A file used for B3 cpm_metadata. Default:
                                --input for single-file mode, first list entry
                                for --input-is-list mode.
  --run-b0-qa                   Run B0 event-index QA. Default: disabled.
  --combined-output PATH        Combined B1/B2/B3 ROOT output. Default:
                                OUT_DIR/PREFIX_B.root
  --no-combined-output          Do not write the combined Job B ROOT file.
  --keep-intermediates          Keep B1/B2 intermediate ROOT files. Default.
  --no-keep-intermediates       Remove B1/B2 after successful B3 and optional
                                combined output.
  --b1-max-pair-dca VALUE       B1 max pair DCA. Default: 2.0
  --b1-min-sin-angle VALUE      B1 minimum sin(opening angle). Default: 1.0e-4
  --b1-max-records VALUE        B1 max records per voxel. Default: 500
  --b1-min-records-per-charge VALUE
                                B1 minimum same-charge records. Default: 2
  --b2-min-entries VALUE        B2 minimum accepted pairs per voxel. Default: 1
  --b2-max-pair-dca VALUE       Optional B2 max pair DCA. Default: -1.0
  --help                        Show this message.

Example:
  scripts/run_cpm_b_chain.sh \
    --input root/Reconstructed/79516/clusters_seeds_79516-0.root_CPMVoxelContainer.root \
    --out-dir root/Reconstructed/79516 \
    --prefix run79516_seg0

  scripts/run_cpm_b_chain.sh \
    --input cpm_filelist.txt --input-is-list \
    --out-dir merged --prefix run79516 \
    --run-b0-qa --no-keep-intermediates
EOF
}

root_string() {
  local value=$1
  value=${value//\\/\\\\}
  value=${value//\"/\\\"}
  printf '"%s"' "$value"
}

first_list_entry() {
  local list_file=$1
  local line
  while IFS= read -r line; do
    if [[ -z "$line" || "${line:0:1}" == "#" ]]; then
      continue
    fi
    printf '%s\n' "$line"
    return 0
  done < "$list_file"
  return 1
}

run_root() {
  local macro_call=$1
  echo
  echo "[run_cpm_b_chain] root -l -b -q ${macro_call}"
  root -l -b -q "$macro_call"
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
MACRO_DIR="${REPO_DIR}/macro"

INPUT=""
INPUT_IS_LIST=0
OUT_DIR="."
PREFIX="CPM"
METADATA=""
B1_MAX_PAIR_DCA="2.0"
B1_MIN_SIN_ANGLE="1.0e-4"
B1_MAX_RECORDS="500"
B1_MIN_RECORDS_PER_CHARGE="2"
B2_MIN_ENTRIES="1"
B2_MAX_PAIR_DCA="-1.0"
RUN_B0_QA=0
WRITE_COMBINED=1
COMBINED_OUTPUT=""
KEEP_INTERMEDIATES=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)
      INPUT=${2:-}
      shift 2
      ;;
    --input-is-list)
      INPUT_IS_LIST=1
      shift
      ;;
    --out-dir)
      OUT_DIR=${2:-}
      shift 2
      ;;
    --prefix)
      PREFIX=${2:-}
      shift 2
      ;;
    --metadata)
      METADATA=${2:-}
      shift 2
      ;;
    --run-b0-qa|--enable-b0-qa)
      RUN_B0_QA=1
      shift
      ;;
    --combined-output)
      COMBINED_OUTPUT=${2:-}
      WRITE_COMBINED=1
      shift 2
      ;;
    --no-combined-output)
      WRITE_COMBINED=0
      shift
      ;;
    --keep-intermediates)
      KEEP_INTERMEDIATES=1
      shift
      ;;
    --no-keep-intermediates)
      KEEP_INTERMEDIATES=0
      shift
      ;;
    --b1-max-pair-dca)
      B1_MAX_PAIR_DCA=${2:-}
      shift 2
      ;;
    --b1-min-sin-angle)
      B1_MIN_SIN_ANGLE=${2:-}
      shift 2
      ;;
    --b1-max-records)
      B1_MAX_RECORDS=${2:-}
      shift 2
      ;;
    --b1-min-records-per-charge)
      B1_MIN_RECORDS_PER_CHARGE=${2:-}
      shift 2
      ;;
    --b2-min-entries)
      B2_MIN_ENTRIES=${2:-}
      shift 2
      ;;
    --b2-max-pair-dca)
      B2_MAX_PAIR_DCA=${2:-}
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      if [[ -z "$INPUT" ]]; then
        INPUT=$1
        shift
      else
        echo "Unknown argument: $1" >&2
        usage >&2
        exit 2
      fi
      ;;
  esac
done

if [[ -z "$INPUT" ]]; then
  echo "Missing --input" >&2
  usage >&2
  exit 2
fi

if [[ ! -e "$INPUT" ]]; then
  echo "Input does not exist: $INPUT" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

if [[ -z "$METADATA" ]]; then
  if [[ "$INPUT_IS_LIST" -eq 1 ]]; then
    METADATA=$(first_list_entry "$INPUT") || {
      echo "Could not find a metadata file in list: $INPUT" >&2
      exit 1
    }
  else
    METADATA=$INPUT
  fi
fi

if [[ ! -e "$METADATA" ]]; then
  echo "Metadata file does not exist: $METADATA" >&2
  exit 1
fi

B0_EVENT_INDEX="${OUT_DIR}/${PREFIX}_B0_event_index.root"
B1_POCA="${OUT_DIR}/${PREFIX}_B1_local_line_poca.root"
B2_CORRECTIONS="${OUT_DIR}/${PREFIX}_B2_voxel_corrections.root"
B3_HISTOGRAMS="${OUT_DIR}/${PREFIX}_B3_average_correction_histograms.root"
if [[ -z "$COMBINED_OUTPUT" ]]; then
  COMBINED_OUTPUT="${OUT_DIR}/${PREFIX}_B.root"
fi
if [[ "$WRITE_COMBINED" -eq 1 ]]; then
  mkdir -p "$(dirname "$COMBINED_OUTPUT")"
fi

INPUT_Q=$(root_string "$INPUT")
B0_Q=$(root_string "$B0_EVENT_INDEX")
B1_Q=$(root_string "$B1_POCA")
B2_Q=$(root_string "$B2_CORRECTIONS")
B3_Q=$(root_string "$B3_HISTOGRAMS")
METADATA_Q=$(root_string "$METADATA")

echo "[run_cpm_b_chain] input: $INPUT"
echo "[run_cpm_b_chain] input_is_list: $INPUT_IS_LIST"
echo "[run_cpm_b_chain] metadata: $METADATA"
echo "[run_cpm_b_chain] output directory: $OUT_DIR"
echo "[run_cpm_b_chain] prefix: $PREFIX"
echo "[run_cpm_b_chain] run_b0_qa: $RUN_B0_QA"
echo "[run_cpm_b_chain] write_combined: $WRITE_COMBINED"
echo "[run_cpm_b_chain] keep_intermediates: $KEEP_INTERMEDIATES"

if [[ "$RUN_B0_QA" -eq 1 ]]; then
  if [[ "$INPUT_IS_LIST" -eq 1 ]]; then
    run_root "${MACRO_DIR}/CPM_B0_BuildEventIndex.C(${INPUT_Q},${B0_Q},true)"
  else
    run_root "${MACRO_DIR}/CPM_B0_BuildEventIndex.C(${INPUT_Q},${B0_Q})"
  fi

  run_root "${MACRO_DIR}/CPM_B0_CheckEventIndex.C(${B0_Q})"
fi

if [[ "$INPUT_IS_LIST" -eq 1 ]]; then
  run_root "${MACRO_DIR}/CPM_B1_LocalLinePoCA.C(${INPUT_Q},${B1_Q},true,${B1_MAX_PAIR_DCA},${B1_MIN_SIN_ANGLE},${B1_MAX_RECORDS},${B1_MIN_RECORDS_PER_CHARGE})"
else
  run_root "${MACRO_DIR}/CPM_B1_LocalLinePoCA.C(${INPUT_Q},${B1_Q},${B1_MAX_PAIR_DCA},${B1_MIN_SIN_ANGLE},${B1_MAX_RECORDS},${B1_MIN_RECORDS_PER_CHARGE})"
fi

run_root "${MACRO_DIR}/CPM_B2_AccumulateVoxelCorrections.C(${B1_Q},${B2_Q},${B2_MIN_ENTRIES},${B2_MAX_PAIR_DCA})"
run_root "${MACRO_DIR}/CPM_B3_WriteAverageCorrectionHistograms.C(${B2_Q},${B3_Q},${METADATA_Q})"
run_root "${MACRO_DIR}/CPM_B3_CheckAverageCorrectionHistograms.C(${B3_Q})"

if [[ "$WRITE_COMBINED" -eq 1 ]]; then
  if ! command -v hadd >/dev/null 2>&1; then
    echo "hadd is required for --combined-output, but it was not found in PATH" >&2
    exit 1
  fi

  combined_inputs=()
  if [[ "$RUN_B0_QA" -eq 1 ]]; then
    combined_inputs+=("$B0_EVENT_INDEX")
  fi
  combined_inputs+=("$B1_POCA" "$B2_CORRECTIONS" "$B3_HISTOGRAMS")

  echo
  echo "[run_cpm_b_chain] hadd -f ${COMBINED_OUTPUT} ${combined_inputs[*]}"
  hadd -f "$COMBINED_OUTPUT" "${combined_inputs[@]}"
fi

if [[ "$KEEP_INTERMEDIATES" -eq 0 ]]; then
  echo
  echo "[run_cpm_b_chain] removing intermediate B1/B2 files"
  rm -f "$B1_POCA" "$B2_CORRECTIONS"
  if [[ -e "$B1_POCA" || -e "$B2_CORRECTIONS" ]]; then
    echo "[run_cpm_b_chain] warning: failed to remove one or more intermediate files" >&2
  fi
fi

echo
echo "[run_cpm_b_chain] done"
if [[ "$RUN_B0_QA" -eq 1 ]]; then
  echo "[run_cpm_b_chain] B0: $B0_EVENT_INDEX"
fi
if [[ "$KEEP_INTERMEDIATES" -eq 1 ]]; then
  echo "[run_cpm_b_chain] B1: $B1_POCA"
  echo "[run_cpm_b_chain] B2: $B2_CORRECTIONS"
else
  echo "[run_cpm_b_chain] B1/B2 intermediates: removed"
fi
echo "[run_cpm_b_chain] B3: $B3_HISTOGRAMS"
if [[ "$WRITE_COMBINED" -eq 1 ]]; then
  echo "[run_cpm_b_chain] combined: $COMBINED_OUTPUT"
fi
