#!/usr/bin/env bash
set -euo pipefail

NEXTFLOW_BIN="${NEXTFLOW_BIN:-/Users/philliprichmondpersonal/.local/bin/nextflow}"
JAVA_CMD="${JAVA_CMD:-/Users/philliprichmondpersonal/.sdkman/candidates/java/current/bin/java}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

TEST_ROOT="${REPO_ROOT}/.test/stub_workflow"
INPUT_DIR="${TEST_ROOT}/input"
WORK_DIR="${TEST_ROOT}/work"
LOG_DIR="${TEST_ROOT}/logs"
TEST_CONFIG="${TEST_ROOT}/test_override.config"

if [[ ! -x "${NEXTFLOW_BIN}" ]]; then
  echo "ERROR: Nextflow executable not found or not executable: ${NEXTFLOW_BIN}" >&2
  exit 1
fi

if [[ -z "${JAVA_CMD}" ]]; then
  if command -v java >/dev/null 2>&1; then
    JAVA_CMD="$(command -v java)"
  elif [[ -n "${JAVA_HOME:-}" && -x "${JAVA_HOME}/bin/java" ]]; then
    JAVA_CMD="${JAVA_HOME}/bin/java"
  fi
fi

if [[ -z "${JAVA_CMD}" || ! -x "${JAVA_CMD}" ]]; then
  echo "ERROR: usable java runtime not found. Set JAVA_CMD=/full/path/to/java or JAVA_HOME." >&2
  exit 127
fi

rm -rf "${TEST_ROOT}"
mkdir -p "${INPUT_DIR}/ubam_dir" "${INPUT_DIR}/aligned" "${INPUT_DIR}/refs" "${INPUT_DIR}/assets" "${WORK_DIR}" "${LOG_DIR}"

UBAM="${INPUT_DIR}/ubam_dir/TEST_STUB.part1.bam"
ALIGNED_BAM="${INPUT_DIR}/aligned/TEST_STUB.t2t.bam"
ALIGNED_BAI="${INPUT_DIR}/aligned/TEST_STUB.t2t.bam.bai"
T2T_FASTA="${INPUT_DIR}/refs/chm13v2.0.fa"
HG38_FASTA="${INPUT_DIR}/refs/hg38_no_alt.fa"
CLINVAR_VCF="${INPUT_DIR}/assets/clinvar.vcf.gz"
CLINVAR_TBI="${INPUT_DIR}/assets/clinvar.vcf.gz.tbi"
SNPEFF_TGZ="${INPUT_DIR}/assets/snpeff_hg38_data.tgz"

touch "${UBAM}" "${ALIGNED_BAM}" "${ALIGNED_BAI}" "${CLINVAR_VCF}" "${CLINVAR_TBI}" "${SNPEFF_TGZ}"
cat > "${T2T_FASTA}" <<'EOF'
>chr4
ACGTACGTACGTACGTACGTACGTACGTACGT
EOF
cat > "${HG38_FASTA}" <<'EOF'
>chr4
ACGTACGTACGTACGTACGTACGTACGTACGT
EOF

cat > "${TEST_CONFIG}" <<'EOF'
docker.enabled = false
conda.enabled = false
wave.enabled = false
process {
  withName: 'MERGE_UNALIGNED_BAMS' { cpus = 1; memory = 1.GB }
  withName: 'ALIGN_READS_TO_T2T' { cpus = 1; memory = 1.GB }
  withName: 'EXTRACT_FSHD_LOCUS' { cpus = 1; memory = 1.GB }
  withName: 'ALIGN_READS_TO_HG38' { cpus = 1; memory = 1.GB }
  withName: 'BLAST_FSHD_READS' { cpus = 1; memory = 1.GB }
  withName: 'PAS_CHECK_FSHD' { cpus = 1; memory = 1.GB }
  withName: 'CLASSIFY_FSHD_READS' { cpus = 1; memory = 1.GB }
  withName: 'EXTRACT_CLASSIFIED_READ_SUBSETS' { cpus = 1; memory = 1.GB }
  withName: 'SUMMARIZE_HAPLOTAGS' { cpus = 1; memory = 1.GB }
  withName: 'PROFILE_FSHD_METHYLATION' { cpus = 1; memory = 1.GB }
  withName: 'CALL_FSHD_VARIANTS' { cpus = 1; memory = 1.GB }
  withName: 'BUILD_FSHD_REPORT' { cpus = 1; memory = 1.GB }
  withName: 'SUMMARIZE_COHORT' { cpus = 1; memory = 1.GB }
}
EOF

run_stub() {
  local mode="$1"
  local out_dir="${TEST_ROOT}/results_${mode}"

  rm -rf "${out_dir}"

  local run_cmd=(
    "${NEXTFLOW_BIN}" run "${REPO_ROOT}/main.nf"
    -stub-run
    -c "${TEST_CONFIG}"
    -work-dir "${WORK_DIR}/${mode}"
    -with-report "${LOG_DIR}/${mode}.report.html"
    -with-trace "${LOG_DIR}/${mode}.trace.txt"
    --sample_id "TEST_${mode^^}"
    --t2t_ref_fasta "${T2T_FASTA}"
    --output_dir "${out_dir}"
  )

  if [[ "${mode}" == "ubam" ]]; then
    run_cmd+=(
      --input_ubam_dir "${INPUT_DIR}/ubam_dir"
    )
  else
    run_cmd+=(
      --input_aligned_bam "${ALIGNED_BAM}"
      --input_aligned_bai "${ALIGNED_BAI}"
      --reuse_input_t2t_alignment false
      --run_variant_calling true
      --hg38_ref_fasta "${HG38_FASTA}"
      --clinvar_vcf_gz "${CLINVAR_VCF}"
      --clinvar_vcf_tbi "${CLINVAR_TBI}"
      --snpeff_data_tgz "${SNPEFF_TGZ}"
    )
  fi

  NXF_HOME="${TEST_ROOT}/.nextflow" \
  NXF_DISABLE_CHECK_LATEST=true \
  JAVA_CMD="${JAVA_CMD}" \
    "${run_cmd[@]}"
}

run_stub "ubam"
run_stub "aligned_variant"

required_outputs=(
  "${TEST_ROOT}/results_ubam/build-fshd-report/TEST_UBAM/TEST_UBAM.fshd.report.html"
  "${TEST_ROOT}/results_ubam/build-fshd-report/TEST_UBAM/TEST_UBAM.fshd.report.summary.tsv"
  "${TEST_ROOT}/results_ubam/summarize-cohort/cohort.fshd.summary.html"
  "${TEST_ROOT}/results_aligned_variant/build-fshd-report/TEST_ALIGNED_VARIANT/TEST_ALIGNED_VARIANT.fshd.report.html"
  "${TEST_ROOT}/results_aligned_variant/call-fshd-variants/TEST_ALIGNED_VARIANT/TEST_ALIGNED_VARIANT.variant.summary.tsv"
  "${TEST_ROOT}/results_aligned_variant/summarize-cohort/cohort.fshd.summary.tsv"
)

for f in "${required_outputs[@]}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: Expected output missing: ${f}" >&2
    exit 1
  fi
done

echo "Workflow stub tests passed."
echo "Outputs:"
echo "  ${TEST_ROOT}/results_ubam"
echo "  ${TEST_ROOT}/results_aligned_variant"
