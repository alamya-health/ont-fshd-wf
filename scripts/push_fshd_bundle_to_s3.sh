#!/usr/bin/env bash
set -euo pipefail

S3_PREFIX_DEFAULT="s3://alamyasingapore-nus-lab-production-processing/reference-genome-and-databases/FSHD_Bundle/"
VARIANT_DIR_DEFAULT="/tmp/ont-fshd-variant-assets"
T2T_DIR_DEFAULT="/tmp/ont-fshd-t2t-assets"
AWS_REGION_DEFAULT="${AWS_REGION:-us-west-1}"

S3_PREFIX="${S3_PREFIX_DEFAULT}"
VARIANT_DIR="${VARIANT_DIR_DEFAULT}"
T2T_DIR="${T2T_DIR_DEFAULT}"
AWS_REGION="${AWS_REGION_DEFAULT}"
DRY_RUN=0

usage() {
  cat <<USAGE
Usage:
  $0 [options]

Uploads the curated FSHD reference bundle assets to S3.

Defaults:
  Variant assets: ${VARIANT_DIR_DEFAULT}
  T2T assets:     ${T2T_DIR_DEFAULT}
  S3 prefix:      ${S3_PREFIX_DEFAULT}

Options:
  --s3-prefix <URI>     Destination S3 prefix
  --variant-dir <DIR>   Local variant asset directory
  --t2t-dir <DIR>       Local T2T asset directory
  --region <REGION>     AWS region for AWS CLI calls (default: ${AWS_REGION_DEFAULT})
  --dry-run             Print uploads without executing them
  -h | --help           Show help
USAGE
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --s3-prefix) S3_PREFIX="${2:-}"; shift 2 ;;
    --variant-dir) VARIANT_DIR="${2:-}"; shift 2 ;;
    --t2t-dir) T2T_DIR="${2:-}"; shift 2 ;;
    --region) AWS_REGION="${2:-}"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage 0 ;;
    *) echo "Unknown option: $1" >&2; usage 1 ;;
  esac
done

command -v aws >/dev/null || { echo "ERROR: aws CLI not found" >&2; exit 1; }

normalize_prefix() {
  local prefix="$1"
  if [[ "${prefix}" != s3://* ]]; then
    echo "ERROR: S3 prefix must start with s3://" >&2
    exit 1
  fi
  if [[ "${prefix}" != */ ]]; then
    prefix="${prefix}/"
  fi
  printf "%s" "${prefix}"
}

require_file() {
  local path="$1"
  if [[ ! -f "${path}" ]]; then
    echo "ERROR: required file not found: ${path}" >&2
    exit 1
  fi
}

S3_PREFIX="$(normalize_prefix "${S3_PREFIX}")"

VARIANT_FILES=(
  "${VARIANT_DIR}/hg38/hg38_no_alt.fa"
  "${VARIANT_DIR}/hg38/hg38_no_alt.fa.fai"
  "${VARIANT_DIR}/clinvar/clinvar.vcf.gz"
  "${VARIANT_DIR}/clinvar/clinvar.vcf.gz.tbi"
  "${VARIANT_DIR}/snpeff/snpeff_hg38_data.tgz"
  "${VARIANT_DIR}/manifest.tsv"
)

T2T_FILES=(
  "${T2T_DIR}/chm13/chm13v2.0.fa"
  "${T2T_DIR}/chm13/chm13v2.0.fa.fai"
  "${T2T_DIR}/chm13/chm13v2.0.mmi"
  "${T2T_DIR}/manifest.tsv"
)

for file in "${VARIANT_FILES[@]}" "${T2T_FILES[@]}"; do
  require_file "${file}"
done

AWS_ARGS=(--region "${AWS_REGION}" --only-show-errors)

upload_file() {
  local src="$1"
  local dest="$2"
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "aws s3 cp ${src} ${dest} ${AWS_ARGS[*]}"
  else
    aws s3 cp "${src}" "${dest}" "${AWS_ARGS[@]}"
  fi
}

echo "Uploading FSHD bundle to ${S3_PREFIX}"

upload_file "${VARIANT_DIR}/hg38/hg38_no_alt.fa" "${S3_PREFIX}variant/hg38/hg38_no_alt.fa"
upload_file "${VARIANT_DIR}/hg38/hg38_no_alt.fa.fai" "${S3_PREFIX}variant/hg38/hg38_no_alt.fa.fai"
if [[ -f "${VARIANT_DIR}/clair3/r1041_e82_400bps_sup_v500.tar.gz" ]]; then
  upload_file "${VARIANT_DIR}/clair3/r1041_e82_400bps_sup_v500.tar.gz" "${S3_PREFIX}variant/clair3/r1041_e82_400bps_sup_v500.tar.gz"
fi
upload_file "${VARIANT_DIR}/clinvar/clinvar.vcf.gz" "${S3_PREFIX}variant/clinvar/clinvar.vcf.gz"
upload_file "${VARIANT_DIR}/clinvar/clinvar.vcf.gz.tbi" "${S3_PREFIX}variant/clinvar/clinvar.vcf.gz.tbi"
upload_file "${VARIANT_DIR}/snpeff/snpeff_hg38_data.tgz" "${S3_PREFIX}variant/snpeff/snpeff_hg38_data.tgz"
upload_file "${VARIANT_DIR}/manifest.tsv" "${S3_PREFIX}variant/manifest.tsv"

upload_file "${T2T_DIR}/chm13/chm13v2.0.fa" "${S3_PREFIX}t2t/chm13/chm13v2.0.fa"
upload_file "${T2T_DIR}/chm13/chm13v2.0.fa.fai" "${S3_PREFIX}t2t/chm13/chm13v2.0.fa.fai"
upload_file "${T2T_DIR}/chm13/chm13v2.0.mmi" "${S3_PREFIX}t2t/chm13/chm13v2.0.mmi"
upload_file "${T2T_DIR}/manifest.tsv" "${S3_PREFIX}t2t/manifest.tsv"

echo "Upload plan complete."
