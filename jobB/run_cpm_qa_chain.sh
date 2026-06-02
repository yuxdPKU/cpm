#!/usr/bin/env bash

if [ -z "${BASH_VERSION:-}" ]; then
  exec /usr/bin/env bash "$0" "$@"
fi

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  jobB/run_cpm_qa_chain.sh --input JOB_A.root [options]
  jobB/run_cpm_qa_chain.sh --input cpm_filelist.txt --input-is-list [options]

Runs the CPM QA/offline diagnostic chain:
  Job A CPMVoxelContainer
  -> QA B0 event index (optional)
  -> QA B1 PoCA pair and batch trees
  -> QA B2 voxel correction tree
  -> QA B3 average-correction histograms
  -> QA B3 histogram check
  -> QA 1D diagnostic plots (default)

Options:
  --input PATH                  Job A CPMVoxelContainer ROOT file, or file list.
  --input-is-list               Treat --input as a newline-separated file list.
  --out-dir DIR                 Output directory. Default: .
  --prefix NAME                 Output filename prefix. Default: CPM_QA
  --metadata PATH               Job A file used for B3 cpm_metadata. Default:
                                --input for single-file mode, first list entry
                                for --input-is-list mode.
  --run-b0-qa                   Also build/check the B0 event index.
  --skip-b0-qa                  Skip B0 event-index QA. Default.
  --max-pair-dca VALUE          B1 max pair DCA. Default: 2.0
  --min-sin-angle VALUE         B1 minimum sin(opening angle) for line solver.
                                Default: 1.0e-4
  --max-records VALUE           B1 max raw records per voxel before skipping.
                                0 disables this safety skip. Default: 0
  --min-records-per-charge VALUE
                                B1 minimum records per charge sign. Default: 2
  --min-pair-pt VALUE           B1 minimum pT for records entering pair loops.
                                Default: 0.5
  --max-pair-records VALUE      B1 max selected records per charge sign in each
                                opposite-charge pair batch. 0 means one
                                unlimited full-voxel batch. Default: 10
  --print-voxel-summaries       Print per-voxel B1 summaries to stdout.
  --no-voxel-summaries          Do not print per-voxel B1 summaries. Default.
  --crossing-solver VALUE       B1 crossing solver: helix or line. Default: helix
  --magnetic-field-z VALUE      Helix-solver Bz field in tesla. Default: 1.4
  --write-pair-tree             Write B1 cpm_poca_pairs.
  --no-pair-tree                Skip the pair-level tree and keep batch summaries. Default.
  --b1-files-per-chunk VALUE    Max Job A files per B1 QA chunk in list mode.
                                0 disables B1 chunking. Default: 100
  --b2-input-mode VALUE         B2 input mode: auto, pairs, or batches. Default: batches
  --b2-max-pair-dca VALUE       Optional B2 pair-row DCA refilter. Default: -1.0
  --min-entries VALUE           Minimum accepted pairs per voxel. Default: 1
  --weighted                    Use pair weights in B2 voxel averaging. Default.
  --unweighted                  Use a simple unweighted average.
  --plots                       Draw QA 1D distributions after the chain. Default.
  --no-plots                    Do not draw QA plots.
  --plot-dir DIR                Plot output directory. Default: OUT_DIR/qa_plots
  --help                        Show this message.

Example:
  jobB/run_cpm_qa_chain.sh \
    --input cpm_filelist.txt --input-is-list \
    --out-dir output/run79516_qa --prefix run79516 --unweighted
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
  echo "[run_cpm_qa_chain] root bool check ${function_call}"
  root -l -b -q -e "gROOT->LoadMacro(${macro_file_q}); bool ok = ${function_call}; gSystem->Exit(ok ? 0 : 1);"
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
MACRO_DIR="${REPO_DIR}/jobB"
PLOT_MACRO="${REPO_DIR}/jobB/plot/qa/CPM_QA_DrawIntermediateDistributions.C"

INPUT=""
INPUT_IS_LIST=0
OUT_DIR="."
PREFIX="CPM_QA"
METADATA=""
RUN_B0_QA=0
B1_MAX_PAIR_DCA="2.0"
B1_MIN_SIN_ANGLE="1.0e-4"
B1_MAX_RECORDS="0"
B1_MIN_RECORDS_PER_CHARGE="2"
B1_PRINT_VOXEL_SUMMARIES=0
B1_MIN_PAIR_PT="0.5"
B1_MAX_PAIR_RECORDS="10"
B1_CROSSING_SOLVER="helix"
B1_MAGNETIC_FIELD_Z="1.4"
B1_WRITE_PAIR_TREE=0
B1_FILES_PER_CHUNK="100"
B2_MIN_ENTRIES="1"
B2_MAX_PAIR_DCA="-1.0"
B2_USE_PAIR_WEIGHTS=1
B2_INPUT_MODE="batches"
B3_MIN_ENTRIES="1"
MAKE_PLOTS=1
PLOT_DIR=""

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
    --run-b0-qa)
      RUN_B0_QA=1
      shift
      ;;
    --skip-b0-qa|--no-b0-qa)
      RUN_B0_QA=0
      shift
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
    --print-voxel-summaries|--b1-print-voxel-summaries)
      B1_PRINT_VOXEL_SUMMARIES=1
      shift
      ;;
    --no-voxel-summaries|--b1-no-voxel-summaries)
      B1_PRINT_VOXEL_SUMMARIES=0
      shift
      ;;
    --crossing-solver|--b1-crossing-solver)
      B1_CROSSING_SOLVER=${2:-}
      shift 2
      ;;
    --magnetic-field-z|--b1-magnetic-field-z)
      B1_MAGNETIC_FIELD_Z=${2:-}
      shift 2
      ;;
    --write-pair-tree)
      B1_WRITE_PAIR_TREE=1
      shift
      ;;
    --no-pair-tree)
      B1_WRITE_PAIR_TREE=0
      shift
      ;;
    --b1-files-per-chunk|--qa-files-per-chunk)
      B1_FILES_PER_CHUNK=${2:-}
      shift 2
      ;;
    --b2-input-mode)
      B2_INPUT_MODE=${2:-}
      shift 2
      ;;
    --b2-max-pair-dca)
      B2_MAX_PAIR_DCA=${2:-}
      shift 2
      ;;
    --min-entries|--b2-min-entries)
      B2_MIN_ENTRIES=${2:-}
      B3_MIN_ENTRIES=${2:-}
      shift 2
      ;;
    --b3-min-entries)
      B3_MIN_ENTRIES=${2:-}
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
    --plots)
      MAKE_PLOTS=1
      shift
      ;;
    --no-plots)
      MAKE_PLOTS=0
      shift
      ;;
    --plot-dir)
      PLOT_DIR=${2:-}
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

if [[ "$B1_CROSSING_SOLVER" != "line" && "$B1_CROSSING_SOLVER" != "helix" ]]; then
  echo "Invalid --crossing-solver: $B1_CROSSING_SOLVER (expected line or helix)" >&2
  exit 2
fi

if [[ "$B2_INPUT_MODE" != "auto" && "$B2_INPUT_MODE" != "pairs" && "$B2_INPUT_MODE" != "batches" ]]; then
  echo "Invalid --b2-input-mode: $B2_INPUT_MODE (expected auto, pairs, or batches)" >&2
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

if [[ -z "$PLOT_DIR" ]]; then
  PLOT_DIR="${OUT_DIR}/qa_plots"
fi

INPUT_Q=$(root_string "$INPUT")
OUT_DIR_Q=$(root_string "$OUT_DIR")
PREFIX_Q=$(root_string "$PREFIX")
B1_CROSSING_SOLVER_Q=$(root_std_string "$B1_CROSSING_SOLVER")
B2_INPUT_MODE_Q=$(root_std_string "$B2_INPUT_MODE")
METADATA_Q=$(root_string "$METADATA")
PLOT_DIR_Q=$(root_string "$PLOT_DIR")

echo "[run_cpm_qa_chain] input: $INPUT"
echo "[run_cpm_qa_chain] input_is_list: $INPUT_IS_LIST"
echo "[run_cpm_qa_chain] metadata: $METADATA"
echo "[run_cpm_qa_chain] output directory: $OUT_DIR"
echo "[run_cpm_qa_chain] prefix: $PREFIX"
echo "[run_cpm_qa_chain] run_b0_qa: $RUN_B0_QA"
echo "[run_cpm_qa_chain] write_pair_tree: $B1_WRITE_PAIR_TREE"
echo "[run_cpm_qa_chain] b1_files_per_chunk: $B1_FILES_PER_CHUNK"
echo "[run_cpm_qa_chain] min_pair_pt: $B1_MIN_PAIR_PT"
echo "[run_cpm_qa_chain] max_pair_records_per_charge_batch: $B1_MAX_PAIR_RECORDS"
echo "[run_cpm_qa_chain] print_voxel_summaries: $B1_PRINT_VOXEL_SUMMARIES"
echo "[run_cpm_qa_chain] crossing_solver: $B1_CROSSING_SOLVER"
echo "[run_cpm_qa_chain] magnetic_field_z: $B1_MAGNETIC_FIELD_Z"
echo "[run_cpm_qa_chain] b2_input_mode: $B2_INPUT_MODE"
echo "[run_cpm_qa_chain] use_pair_weights: $B2_USE_PAIR_WEIGHTS"
echo "[run_cpm_qa_chain] min_entries_per_voxel: $B2_MIN_ENTRIES"
echo "[run_cpm_qa_chain] make_plots: $MAKE_PLOTS"

run_root_bool_check "${MACRO_DIR}/CPM_QA_RunOfflineDiagnostics.C" "CPM_QA_RunOfflineDiagnostics(${INPUT_Q},${OUT_DIR_Q},${PREFIX_Q},${INPUT_IS_LIST},${RUN_B0_QA},${B1_MAX_PAIR_DCA},${B1_MIN_SIN_ANGLE},${B1_MAX_RECORDS},${B1_MIN_RECORDS_PER_CHARGE},${B1_PRINT_VOXEL_SUMMARIES},${B1_MIN_PAIR_PT},${B1_MAX_PAIR_RECORDS},${B1_CROSSING_SOLVER_Q},${B1_MAGNETIC_FIELD_Z},${B1_WRITE_PAIR_TREE},${B2_MIN_ENTRIES},${B2_MAX_PAIR_DCA},${B2_USE_PAIR_WEIGHTS},${B2_INPUT_MODE_Q},${B3_MIN_ENTRIES},${METADATA_Q},${B1_FILES_PER_CHUNK})"

if [[ "$MAKE_PLOTS" -eq 1 ]]; then
  run_root_bool_check "$PLOT_MACRO" "CPM_QA_DrawIntermediateDistributions(${OUT_DIR_Q},${PREFIX_Q},${PLOT_DIR_Q},true)"
fi

echo
echo "[run_cpm_qa_chain] done"
echo "[run_cpm_qa_chain] B1: ${OUT_DIR}/${PREFIX}_QA_B1_poca.root or ${OUT_DIR}/${PREFIX}_QA_B1_chunks/"
echo "[run_cpm_qa_chain] B2: ${OUT_DIR}/${PREFIX}_QA_B2_voxel_corrections.root"
echo "[run_cpm_qa_chain] B3: ${OUT_DIR}/${PREFIX}_B3_average_correction_histograms.root"
if [[ "$RUN_B0_QA" -eq 1 ]]; then
  echo "[run_cpm_qa_chain] B0: ${OUT_DIR}/${PREFIX}_QA_B0_event_index.root"
fi
if [[ "$MAKE_PLOTS" -eq 1 ]]; then
  echo "[run_cpm_qa_chain] plots: $PLOT_DIR"
fi
