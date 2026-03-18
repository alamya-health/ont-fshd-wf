#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-/tmp/ont-fshd-t2t-assets}"
THREADS="${THREADS:-4}"
HS1_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"
SAMTOOLS_IMAGE="quay.io/biocontainers/samtools:1.22.1--h96c455f_0"
MINIMAP2_IMAGE="quay.io/biocontainers/minimap2:2.28--he4a0461_3"
MINIMAP2_BIN="${MINIMAP2_BIN:-}"

mkdir -p "${OUT_DIR}"/{downloads,chm13}

fetch() {
  local url="$1"
  local out="$2"
  if [[ -s "${out}" ]]; then
    echo "Already exists: ${out}"
    return 0
  fi
  curl -L --fail --silent "${url}" -o "${out}"
}

echo "Downloading T2T-CHM13v2.0 FASTA..."
fetch "${HS1_URL}" "${OUT_DIR}/downloads/hs1.fa.gz"

if [[ ! -s "${OUT_DIR}/chm13/chm13v2.0.fa" ]]; then
  echo "Expanding T2T-CHM13v2.0 FASTA..."
  gzip -dc "${OUT_DIR}/downloads/hs1.fa.gz" > "${OUT_DIR}/chm13/chm13v2.0.fa"
fi

if [[ ! -s "${OUT_DIR}/chm13/chm13v2.0.fa.fai" ]]; then
  echo "Indexing T2T-CHM13v2.0 FASTA with samtools..."
  docker run --rm -v "${OUT_DIR}/chm13:/data" "${SAMTOOLS_IMAGE}" \
    samtools faidx /data/chm13v2.0.fa
fi

if [[ ! -s "${OUT_DIR}/chm13/chm13v2.0.mmi" ]]; then
  echo "Building minimap2 index for T2T-CHM13v2.0..."
  if [[ -n "${MINIMAP2_BIN}" ]]; then
    "${MINIMAP2_BIN}" -t "${THREADS}" -d "${OUT_DIR}/chm13/chm13v2.0.mmi" "${OUT_DIR}/chm13/chm13v2.0.fa"
  elif command -v minimap2 >/dev/null 2>&1; then
    minimap2 -t "${THREADS}" -d "${OUT_DIR}/chm13/chm13v2.0.mmi" "${OUT_DIR}/chm13/chm13v2.0.fa"
  else
    docker run --rm -v "${OUT_DIR}/chm13:/data" "${MINIMAP2_IMAGE}" \
      minimap2 -t "${THREADS}" -d /data/chm13v2.0.mmi /data/chm13v2.0.fa
  fi
fi

MANIFEST="${OUT_DIR}/manifest.tsv"
{
  printf "asset_name\tpath\tsha256\tsize_bytes\n"
  for file in \
    "${OUT_DIR}/chm13/chm13v2.0.fa" \
    "${OUT_DIR}/chm13/chm13v2.0.fa.fai" \
    "${OUT_DIR}/chm13/chm13v2.0.mmi"
  do
    printf "%s\t%s\t%s\t%s\n" \
      "$(basename "${file}")" \
      "${file}" \
      "$(shasum -a 256 "${file}" | awk '{print $1}')" \
      "$(stat -f%z "${file}")"
  done
} > "${MANIFEST}"

echo "Asset staging complete:"
echo "  ${OUT_DIR}"
echo "Manifest:"
echo "  ${MANIFEST}"
