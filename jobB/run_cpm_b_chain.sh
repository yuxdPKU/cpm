#!/usr/bin/env bash

if [ -z "${BASH_VERSION:-}" ]; then
  exec /usr/bin/env bash "$0" "$@"
fi

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  jobB/run_cpm_b_chain.sh --input JOB_A.root [options]
  jobB/run_cpm_b_chain.sh --input cpm_filelist.txt --input-is-list [options]

Runs the CPM production Job B calculation:
  Job A CPMVoxelContainer
  -> opposite-charge crossing-point PoCA
  -> voxel averaging
  -> B3 average-correction histogram output
  B3 histogram check

For pair trees, batch trees, voxel summaries, or B0 event-index diagnostics,
run jobB/CPM_QA_RunOfflineDiagnostics.C instead.

Options:
  --input PATH                  Job A CPMVoxelContainer ROOT file, or file list.
  --input-is-list               Treat --input as a newline-separated file list.
  --out-dir DIR                 Output directory. Default: .
  --prefix NAME                 Output filename prefix. Default: CPM
  --metadata PATH               Job A file used for B3 cpm_metadata. Default:
                                --input for single-file mode, first list entry
                                for --input-is-list mode.
  --max-pair-dca VALUE          Max pair DCA. Default: 2.0
  --min-sin-angle VALUE         Minimum sin(opening angle) for line solver.
                                Default: 1.0e-4
  --max-records VALUE           Max raw records per voxel before skipping.
                                0 disables this safety skip. Default: 0
  --min-records-per-charge VALUE
                                Minimum records per charge sign. Default: 2
  --min-pair-pt VALUE           Minimum pT for records entering pair loops.
                                Default: 0.5
  --max-pair-records VALUE      Max selected records per charge sign in each
                                opposite-charge pair batch after keeping the
                                closest record per unique track. 0 means one
                                unlimited full-voxel batch.
                                Default: 10
  --crossing-solver VALUE        Crossing solver: helix or line. Default: helix
  --magnetic-field-z VALUE       Helix-solver Bz field in tesla. Default: 1.4
  --min-entries VALUE           Minimum accepted pairs per voxel. Default: 1
  --weighted                    Use pair weights in voxel averaging. Default.
  --unweighted                  Use a simple unweighted average.
  --help                        Show this message.

Example:
  jobB/run_cpm_b_chain.sh \
    --input root/Reconstructed/79516/clusters_seeds_79516-0.root_CPMVoxelContainer.root \
    --out-dir root/Reconstructed/79516 \
    --prefix run79516_seg0

  jobB/run_cpm_b_chain.sh \
    --input cpm_filelist.txt --input-is-list \
    --out-dir merged --prefix run79516 \
    --unweighted
EOF
}

root_string() {
  local value=$1
  value=${value//\\/\\\\}
  value=${value//\"/\\\"}
  printf '"%s"' "$value"
}

root_std_string() {
  printf 'std::string(%s)' "$(root_string "$1")"
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

run_root_bool_check() {
  local macro_file=$1
  local function_call=$2
  local macro_file_q
  macro_file_q=$(root_string "$macro_file")
  echo
  echo "[run_cpm_b_chain] root bool check ${function_call}"
  root -l -b -q -e "gROOT->LoadMacro(${macro_file_q}); bool ok = ${function_call}; gSystem->Exit(ok ? 0 : 1);"
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
MACRO_DIR="${REPO_DIR}/jobB"

INPUT=""
INPUT_IS_LIST=0
OUT_DIR="."
PREFIX="CPM"
METADATA=""
B1_MAX_PAIR_DCA="2.0"
B1_MIN_SIN_ANGLE="1.0e-4"
B1_MAX_RECORDS="0"
B1_MIN_RECORDS_PER_CHARGE="2"
B1_MIN_PAIR_PT="0.5"
B1_MAX_PAIR_RECORDS="10"
B1_CROSSING_SOLVER="helix"
B1_MAGNETIC_FIELD_Z="1.4"
B2_MIN_ENTRIES="1"
B2_USE_PAIR_WEIGHTS=1

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
    --max-pair-dca|--b1-max-pair-dca)
      B1_MAX_PAIR_DCA=${2:-}
      shift 2
      ;;
    --min-sin-angle|--b1-min-sin-angle)
      B1_MIN_SIN_ANGLE=${2:-}
      shift 2
      ;;
    --max-records|--b1-max-records)
      B1_MAX_RECORDS=${2:-}
      shift 2
      ;;
    --min-records-per-charge|--b1-min-records-per-charge)
      B1_MIN_RECORDS_PER_CHARGE=${2:-}
      shift 2
      ;;
    --min-pair-pt|--b1-min-pair-pt)
      B1_MIN_PAIR_PT=${2:-}
      shift 2
      ;;
    --max-pair-records|--b1-max-pair-records)
      B1_MAX_PAIR_RECORDS=${2:-}
      shift 2
      ;;
    --crossing-solver|--b1-crossing-solver)
      B1_CROSSING_SOLVER=${2:-}
      shift 2
      ;;
    --magnetic-field-z|--b1-magnetic-field-z)
      B1_MAGNETIC_FIELD_Z=${2:-}
      shift 2
      ;;
    --min-entries|--b2-min-entries)
      B2_MIN_ENTRIES=${2:-}
      shift 2
      ;;
    --weighted|--b2-weighted)
      B2_USE_PAIR_WEIGHTS=1
      shift
      ;;
    --unweighted|--b2-unweighted)
      B2_USE_PAIR_WEIGHTS=0
      shift
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

if [[ "$B1_CROSSING_SOLVER" != "line" && "$B1_CROSSING_SOLVER" != "helix" ]]; then
  echo "Invalid --crossing-solver: $B1_CROSSING_SOLVER (expected line or helix)" >&2
  exit 2
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

B3_HISTOGRAMS="${OUT_DIR}/${PREFIX}_B3_average_correction_histograms.root"

INPUT_Q=$(root_string "$INPUT")
B3_Q=$(root_string "$B3_HISTOGRAMS")
METADATA_Q=$(root_string "$METADATA")
B1_CROSSING_SOLVER_Q=$(root_std_string "$B1_CROSSING_SOLVER")

echo "[run_cpm_b_chain] input: $INPUT"
echo "[run_cpm_b_chain] input_is_list: $INPUT_IS_LIST"
echo "[run_cpm_b_chain] metadata: $METADATA"
echo "[run_cpm_b_chain] output directory: $OUT_DIR"
echo "[run_cpm_b_chain] prefix: $PREFIX"
echo "[run_cpm_b_chain] min_pair_pt: $B1_MIN_PAIR_PT"
echo "[run_cpm_b_chain] max_pair_records_per_charge_batch: $B1_MAX_PAIR_RECORDS"
echo "[run_cpm_b_chain] crossing_solver: $B1_CROSSING_SOLVER"
echo "[run_cpm_b_chain] magnetic_field_z: $B1_MAGNETIC_FIELD_Z"
echo "[run_cpm_b_chain] use_pair_weights: $B2_USE_PAIR_WEIGHTS"
echo "[run_cpm_b_chain] min_entries_per_voxel: $B2_MIN_ENTRIES"

run_root_bool_check "${MACRO_DIR}/CPM_ComputeAverageCorrection.C" "CPM_ComputeAverageCorrection(${INPUT_Q},${B3_Q},${INPUT_IS_LIST},${B2_USE_PAIR_WEIGHTS},${B2_MIN_ENTRIES},${B1_MAX_PAIR_DCA},${B1_MIN_SIN_ANGLE},${B1_MAX_RECORDS},${B1_MIN_RECORDS_PER_CHARGE},${B1_MIN_PAIR_PT},${B1_MAX_PAIR_RECORDS},${B1_CROSSING_SOLVER_Q},${B1_MAGNETIC_FIELD_Z},${METADATA_Q})"

run_root_bool_check "${MACRO_DIR}/CPM_B3_CheckAverageCorrectionHistograms.C" "CPM_B3_CheckAverageCorrectionHistograms(${B3_Q})"

echo
echo "[run_cpm_b_chain] done"
echo "[run_cpm_b_chain] B3: $B3_HISTOGRAMS"
