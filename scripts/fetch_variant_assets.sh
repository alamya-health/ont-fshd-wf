#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-/tmp/ont-fshd-variant-assets}"
CLAIR3_MODEL="${CLAIR3_MODEL:-r1041_e82_400bps_sup_v500}"
HG38_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
CLAIR3_URL="https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/${CLAIR3_MODEL}.tar.gz"
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
SNPEFF_URL="https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip"

mkdir -p "${OUT_DIR}"/{downloads,hg38,clinvar,clair3,snpeff,tmp}

fetch() {
  local url="$1"
  local out="$2"
  if [[ -s "${out}" ]]; then
    echo "Already exists: ${out}"
    return 0
  fi
  curl -L --fail --silent "${url}" -o "${out}"
}

echo "Downloading HG38 no-alt FASTA..."
fetch "${HG38_URL}" "${OUT_DIR}/downloads/hg38_no_alt.fa.gz"

echo "Downloading Clair3 model..."
fetch "${CLAIR3_URL}" "${OUT_DIR}/clair3/${CLAIR3_MODEL}.tar.gz"

echo "Downloading ClinVar..."
fetch "${CLINVAR_URL}" "${OUT_DIR}/clinvar/clinvar.vcf.gz"
fetch "${CLINVAR_URL}.tbi" "${OUT_DIR}/clinvar/clinvar.vcf.gz.tbi"

echo "Downloading snpEff core..."
fetch "${SNPEFF_URL}" "${OUT_DIR}/downloads/snpEff_latest_core.zip"

if [[ ! -s "${OUT_DIR}/hg38/hg38_no_alt.fa" ]]; then
  echo "Expanding HG38 FASTA..."
  gzip -dc "${OUT_DIR}/downloads/hg38_no_alt.fa.gz" > "${OUT_DIR}/hg38/hg38_no_alt.fa"
fi

if [[ ! -s "${OUT_DIR}/hg38/hg38_no_alt.fa.fai" ]]; then
  echo "Indexing HG38 FASTA with samtools..."
  docker run --rm -v "${OUT_DIR}/hg38:/data" quay.io/biocontainers/samtools:1.22.1--h96c455f_0 \
    samtools faidx /data/hg38_no_alt.fa
fi

if [[ ! -s "${OUT_DIR}/snpeff/snpEff_v_latest_core/snpEff.jar" ]]; then
  echo "Expanding snpEff core..."
  unzip -q -o "${OUT_DIR}/downloads/snpEff_latest_core.zip" -d "${OUT_DIR}/snpeff"
fi

SNPEFF_HOME=""
if [[ -f "${OUT_DIR}/snpeff/snpEff/snpEff.jar" ]]; then
  SNPEFF_HOME="${OUT_DIR}/snpeff/snpEff"
elif [[ -f "${OUT_DIR}/snpeff/snpEff_v_latest_core/snpEff.jar" ]]; then
  SNPEFF_HOME="${OUT_DIR}/snpeff/snpEff_v_latest_core"
else
  SNPEFF_HOME="$(find "${OUT_DIR}/snpeff" -maxdepth 3 -type f -name snpEff.jar -print | head -n 1 | xargs dirname)"
fi

if [[ -z "${SNPEFF_HOME}" || ! -f "${SNPEFF_HOME}/snpEff.jar" ]]; then
  echo "Could not locate snpEff.jar under ${OUT_DIR}/snpeff" >&2
  exit 1
fi

if [[ ! -d "${SNPEFF_HOME}/data/hg38" && ! -d "${SNPEFF_HOME}/hg38" ]]; then
  echo "Downloading snpEff hg38 database..."
  docker run --rm \
    -v "${SNPEFF_HOME}:/opt/snpeff" \
    eclipse-temurin:21-jre \
    java -jar /opt/snpeff/snpEff.jar download -dataDir /opt/snpeff hg38
fi

if [[ ! -s "${OUT_DIR}/snpeff/snpeff_hg38_data.tgz" ]] || ! tar -tzf "${OUT_DIR}/snpeff/snpeff_hg38_data.tgz" 2>/dev/null | grep -q '^data/hg38/'; then
  echo "Packaging snpEff hg38 data..."
  STAGE_DIR="${OUT_DIR}/tmp/snpeff-data-stage"
  rm -rf "${STAGE_DIR}"
  mkdir -p "${STAGE_DIR}/data"
  if [[ -d "${SNPEFF_HOME}/data/hg38" ]]; then
    cp -R "${SNPEFF_HOME}/data/hg38" "${STAGE_DIR}/data/hg38"
  elif [[ -d "${SNPEFF_HOME}/hg38" ]]; then
    cp -R "${SNPEFF_HOME}/hg38" "${STAGE_DIR}/data/hg38"
  else
    echo "Could not locate snpEff hg38 database under ${SNPEFF_HOME}" >&2
    exit 1
  fi
  tar -C "${STAGE_DIR}" -czf "${OUT_DIR}/snpeff/snpeff_hg38_data.tgz" data
  rm -rf "${STAGE_DIR}"
fi

MANIFEST="${OUT_DIR}/manifest.tsv"
{
  printf "asset_name\tpath\tsha256\tsize_bytes\n"
  for file in \
    "${OUT_DIR}/hg38/hg38_no_alt.fa" \
    "${OUT_DIR}/hg38/hg38_no_alt.fa.fai" \
    "${OUT_DIR}/clair3/${CLAIR3_MODEL}.tar.gz" \
    "${OUT_DIR}/clinvar/clinvar.vcf.gz" \
    "${OUT_DIR}/clinvar/clinvar.vcf.gz.tbi" \
    "${OUT_DIR}/snpeff/snpeff_hg38_data.tgz"
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
