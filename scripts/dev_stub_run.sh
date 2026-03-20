#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEV_DIR="${ROOT_DIR}/.dev/stub_run"
INPUT_DIR="${DEV_DIR}/inputs"
RESULTS_DIR="${DEV_DIR}/results"
OVERRIDE_CONFIG="${DEV_DIR}/stub_override.config"
NXF_HOME_DIR="${DEV_DIR}/.nextflow"
NEXTFLOW_BIN="${NEXTFLOW_BIN:-/Users/philliprichmondpersonal/.local/bin/nextflow}"
JAVA_CMD="${JAVA_CMD:-/Users/philliprichmondpersonal/.sdkman/candidates/java/current/bin/java}"

mkdir -p "${INPUT_DIR}/ubam_dir" "${RESULTS_DIR}" "${NXF_HOME_DIR}"

if [[ ! -x "${NEXTFLOW_BIN}" ]]; then
  echo "ERROR: nextflow not found or not executable at ${NEXTFLOW_BIN}"
  exit 127
fi

if [[ -z "${JAVA_CMD}" ]]; then
  if command -v java >/dev/null 2>&1; then
    JAVA_CMD="$(command -v java)"
  elif [[ -n "${JAVA_HOME:-}" && -x "${JAVA_HOME}/bin/java" ]]; then
    JAVA_CMD="${JAVA_HOME}/bin/java"
  fi
fi

if [[ -z "${JAVA_CMD}" || ! -x "${JAVA_CMD}" ]]; then
  echo "ERROR: usable java runtime not found. Set JAVA_CMD=/full/path/to/java or JAVA_HOME."
  exit 127
fi

UBAM="${INPUT_DIR}/ubam_dir/DEV_STUB.part1.bam"
T2T_FASTA="${INPUT_DIR}/chm13v2.0.fa"

touch "${UBAM}"
cat > "${T2T_FASTA}" <<'EOF'
>chr4
ACGTACGTACGTACGTACGTACGTACGTACGT
EOF

cat > "${OVERRIDE_CONFIG}" <<'EOF'
docker.enabled = false
conda.enabled = false
wave.enabled = false
EOF

NXF_HOME="${NXF_HOME_DIR}" \
NXF_DISABLE_CHECK_LATEST=true \
JAVA_CMD="${JAVA_CMD}" "${NEXTFLOW_BIN}" run "${ROOT_DIR}/main.nf" \
  -stub-run \
  -resume \
  -c "${OVERRIDE_CONFIG}" \
  --sample_id DEV_STUB \
  --input_ubam_dir "${INPUT_DIR}/ubam_dir" \
  --t2t_ref_fasta "${T2T_FASTA}" \
  --output_dir "${RESULTS_DIR}"
