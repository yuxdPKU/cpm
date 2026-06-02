#!/usr/bin/env bash

if [ -z "${BASH_VERSION:-}" ]; then
  exec /usr/bin/env bash "$0" "$@"
fi

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  jobB/run_cpm_closure_slice.sh --qa-dir QA_DIR --prefix PREFIX [options]

Draws fixed-phi/fixed-z CPM closure diagnostics from QA B1 cpm_poca_pairs.
The output directory is separate from the regular qa_plots directory.

Options:
  --qa-dir DIR             QA output directory containing PREFIX_QA_B1_chunks.txt
                           or PREFIX_QA_B1_poca.root.
  --prefix NAME            QA filename prefix.
  --plot-dir DIR           Output plot directory. Default:
                           QA_DIR/closure_slice_phiPHI_zZ
  --phi VALUE              Selected phi in radians. Default: 4.7
  --z VALUE                Selected |z| in cm. Default: 10.0
  --dca-thresholds CSV     DCA thresholds to refilter pair rows. Default:
                           2.0,1.0,0.5,0.2
  --b1-input PATH          Explicit B1 ROOT file or B1 chunk-list file.
  --b1-input-is-list       Treat --b1-input as a list of B1 ROOT files.
  --b3-file PATH           Optional B3 histogram file to overlay. Default:
                           QA_DIR/PREFIX_B3_average_correction_histograms.root
  --help                   Show this message.

Example:
  jobB/run_cpm_closure_slice.sh \
    --qa-dir jobB/output/sim_genfit_unweighted_qa \
    --prefix sim_genfit_unweighted \
    --phi 4.7 --z 10 \
    --dca-thresholds 2.0,1.0,0.5,0.2
EOF
}

root_string() {
  local value=$1
  value=${value//\\/\\\\}
  value=${value//\"/\\\"}
  printf '"%s"' "$value"
}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
MACRO="${REPO_DIR}/jobB/plot/qa/CPM_QA_DrawClosureSlice.C"

QA_DIR=""
PREFIX=""
PLOT_DIR=""
PHI="4.7"
SELECT_Z="10.0"
DCA_THRESHOLDS="2.0,1.0,0.5,0.2"
B1_INPUT=""
B1_INPUT_IS_LIST=0
B3_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --qa-dir)
      QA_DIR=${2:-}
      shift 2
      ;;
    --prefix)
      PREFIX=${2:-}
      shift 2
      ;;
    --plot-dir)
      PLOT_DIR=${2:-}
      shift 2
      ;;
    --phi)
      PHI=${2:-}
      shift 2
      ;;
    --z|--select-z)
      SELECT_Z=${2:-}
      shift 2
      ;;
    --dca-thresholds)
      DCA_THRESHOLDS=${2:-}
      shift 2
      ;;
    --b1-input)
      B1_INPUT=${2:-}
      shift 2
      ;;
    --b1-input-is-list)
      B1_INPUT_IS_LIST=1
      shift
      ;;
    --b3-file)
      B3_FILE=${2:-}
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$QA_DIR" || -z "$PREFIX" ]]; then
  echo "Missing --qa-dir or --prefix" >&2
  usage >&2
  exit 2
fi

if [[ ! -e "$MACRO" ]]; then
  echo "Missing macro: $MACRO" >&2
  exit 1
fi

if [[ -z "$B1_INPUT" && ! -d "$QA_DIR" ]]; then
  echo "QA directory does not exist: $QA_DIR" >&2
  exit 1
fi

if [[ -n "$B1_INPUT" && ! -e "$B1_INPUT" ]]; then
  echo "B1 input does not exist: $B1_INPUT" >&2
  exit 1
fi

MACRO_Q=$(root_string "$MACRO")
QA_DIR_Q=$(root_string "$QA_DIR")
PREFIX_Q=$(root_string "$PREFIX")
PLOT_DIR_Q=$(root_string "$PLOT_DIR")
DCA_Q=$(root_string "$DCA_THRESHOLDS")
B1_INPUT_Q=$(root_string "$B1_INPUT")
B3_FILE_Q=$(root_string "$B3_FILE")

FUNCTION_CALL="CPM_QA_DrawClosureSlice(${QA_DIR_Q},${PREFIX_Q},${PLOT_DIR_Q},${PHI},${SELECT_Z},${DCA_Q},${B1_INPUT_Q},${B1_INPUT_IS_LIST},${B3_FILE_Q})"

echo "[run_cpm_closure_slice] qa_dir: $QA_DIR"
echo "[run_cpm_closure_slice] prefix: $PREFIX"
echo "[run_cpm_closure_slice] plot_dir: ${PLOT_DIR:-<macro default>}"
echo "[run_cpm_closure_slice] phi: $PHI"
echo "[run_cpm_closure_slice] |z|: $SELECT_Z"
echo "[run_cpm_closure_slice] dca thresholds: $DCA_THRESHOLDS"
if [[ -n "$B1_INPUT" ]]; then
  echo "[run_cpm_closure_slice] b1_input: $B1_INPUT"
  echo "[run_cpm_closure_slice] b1_input_is_list: $B1_INPUT_IS_LIST"
fi
if [[ -n "$B3_FILE" ]]; then
  echo "[run_cpm_closure_slice] b3_file: $B3_FILE"
fi

echo
echo "[run_cpm_closure_slice] root bool check ${FUNCTION_CALL}"
root -l -b -q -e "gROOT->LoadMacro(${MACRO_Q}); bool ok = ${FUNCTION_CALL}; gSystem->Exit(ok ? 0 : 1);"
