#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Render the FSHD HTML report from already-published workflow outputs.

Usage:
  scripts/render_report_from_outputs.sh \
    --sample-id SAMPLE \
    --output-root /path/to/published/output/root \
    [--outdir /path/to/write/report] \
    [--contraction-threshold 8] \
    [--locus-bed assets/fshd_locus_t2t.bed] \
    [--methyl-bed assets/fshd_methylation_regions_t2t.bed]

Expected layout under --output-root:
  classify-fshd-reads/<sample>/<sample>.fshd.classification/
  extract-classified-read-subsets/<sample>/<sample>.classified-subsets/
  extract-fshd-locus/<sample>/<sample>.fshd.locus.flagstat.txt
  extract-fshd-locus/<sample>/<sample>.fshd.locus.coverage.tsv
  summarize-haplotags/<sample>/<sample>.haplotag.summary.tsv
  profile-fshd-methylation/<sample>/<sample>.methylation/
  profile-fshd-methylation/<sample>/<sample>.methylation.summary.tsv

The script writes:
  <outdir>/<sample>.fshd.report.html
  <outdir>/<sample>.fshd.report.summary.tsv

It also accepts a flat scratch directory containing copied artifacts like:
  A_4qA_all-reads.csv
  4qA_reads_blast.csv
  <sample>.classified-subsets.manifest.tsv
  <sample>.methylation.summary.tsv
  4qA_all.methyl.bed.gz
EOF
}

sample_id=""
output_root=""
outdir=""
contraction_threshold="8"
locus_bed="assets/fshd_locus_t2t.bed"
methyl_bed="assets/fshd_methylation_regions_t2t.bed"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) sample_id="${2:-}"; shift 2 ;;
    --output-root) output_root="${2:-}"; shift 2 ;;
    --outdir) outdir="${2:-}"; shift 2 ;;
    --contraction-threshold) contraction_threshold="${2:-}"; shift 2 ;;
    --locus-bed) locus_bed="${2:-}"; shift 2 ;;
    --methyl-bed) methyl_bed="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$sample_id" || -z "$output_root" ]]; then
  usage
  exit 1
fi

if [[ -z "$outdir" ]]; then
  outdir="${output_root}/report-refresh/${sample_id}"
fi

classification_dir="${output_root}/classify-fshd-reads/${sample_id}/${sample_id}.fshd.classification"
subset_dir="${output_root}/extract-classified-read-subsets/${sample_id}/${sample_id}.classified-subsets"
flagstat_txt="${output_root}/extract-fshd-locus/${sample_id}/${sample_id}.fshd.locus.flagstat.txt"
coverage_tsv="${output_root}/extract-fshd-locus/${sample_id}/${sample_id}.fshd.locus.coverage.tsv"
hap_tsv="${output_root}/summarize-haplotags/${sample_id}/${sample_id}.haplotag.summary.tsv"
methylation_dir="${output_root}/profile-fshd-methylation/${sample_id}/${sample_id}.methylation"
methyl_tsv="${output_root}/profile-fshd-methylation/${sample_id}/${sample_id}.methylation.summary.tsv"

if [[ ! -d "$classification_dir" && -f "${output_root}/A_4qA_all-reads.csv" ]]; then
  tmp_root="$(mktemp -d /tmp/fshd-report-flat.XXXXXX)"
  classification_dir="${tmp_root}/classification"
  subset_dir="${tmp_root}/subsets"
  methylation_dir="${tmp_root}/methylation"
  mkdir -p "${classification_dir}/blast_results" "$subset_dir" "$methylation_dir"

  cp "${output_root}/A_4qA_all-reads.csv" "${classification_dir}/"
  [[ -f "${output_root}/A_4qA_complete-reads.csv" ]] && cp "${output_root}/A_4qA_complete-reads.csv" "${classification_dir}/"
  [[ -f "${output_root}/B_4qB_all-reads.csv" ]] && cp "${output_root}/B_4qB_all-reads.csv" "${classification_dir}/"
  [[ -f "${output_root}/B_4qB_complete-reads.csv" ]] && cp "${output_root}/B_4qB_complete-reads.csv" "${classification_dir}/"
  [[ -f "${output_root}/chr10_all-reads.csv" ]] && cp "${output_root}/chr10_all-reads.csv" "${classification_dir}/"
  [[ -f "${output_root}/chr10_complete-reads.csv" ]] && cp "${output_root}/chr10_complete-reads.csv" "${classification_dir}/"
  [[ -f "${output_root}/chimeric_reads.csv" ]] && cp "${output_root}/chimeric_reads.csv" "${classification_dir}/"
  [[ -f "${output_root}/chr4_undefined_all-reads.csv" ]] && cp "${output_root}/chr4_undefined_all-reads.csv" "${classification_dir}/"
  [[ -f "${output_root}/4qA_reads_blast.csv" ]] && cp "${output_root}/4qA_reads_blast.csv" "${classification_dir}/blast_results/"

  if [[ -f "${output_root}/${sample_id}.classified-subsets.manifest.tsv" ]]; then
    cp "${output_root}/${sample_id}.classified-subsets.manifest.tsv" "${subset_dir}/"
  fi

  if [[ -f "${output_root}/${sample_id}.methylation.summary.tsv" ]]; then
    methyl_tsv="${output_root}/${sample_id}.methylation.summary.tsv"
  fi
  [[ -f "${output_root}/4qA_all.methyl.bed.gz" ]] && cp "${output_root}/4qA_all.methyl.bed.gz" "${methylation_dir}/"
  [[ -f "${output_root}/4qA_all.stats.tsv" ]] && cp "${output_root}/4qA_all.stats.tsv" "${methylation_dir}/"

  if [[ -f "${output_root}/${sample_id}.fshd.report.html" ]]; then
    :
  fi

  if [[ ! -f "${output_root}/flagstat.txt" ]]; then
    cat <<EOF > "${tmp_root}/flagstat.txt"
flagstat unavailable in flat input bundle
EOF
    flagstat_txt="${tmp_root}/flagstat.txt"
  fi

  if [[ ! -f "${output_root}/coverage.tsv" ]]; then
    cat <<EOF > "${tmp_root}/coverage.tsv"
#chrom	start	end	label	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
EOF
    coverage_tsv="${tmp_root}/coverage.tsv"
  fi

  if [[ ! -f "${output_root}/haplotag.summary.tsv" && ! -f "$hap_tsv" ]]; then
    cat <<EOF > "${tmp_root}/hap.tsv"
subset_name	hp	count	total_subset_reads	has_hp_tags
4qA_all	UNSET	0	0	false
EOF
    hap_tsv="${tmp_root}/hap.tsv"
  elif [[ -f "${output_root}/haplotag.summary.tsv" ]]; then
    hap_tsv="${output_root}/haplotag.summary.tsv"
  fi

fi

required_paths=(
  "$classification_dir"
  "$subset_dir"
  "$flagstat_txt"
  "$coverage_tsv"
  "$hap_tsv"
  "$methylation_dir"
  "$methyl_tsv"
  "$locus_bed"
  "$methyl_bed"
)

for path in "${required_paths[@]}"; do
  if [[ ! -e "$path" ]]; then
    echo "Missing required report input: $path" >&2
    exit 1
  fi
done

mkdir -p "$outdir"

PYTHONDONTWRITEBYTECODE=1 python3 "${PWD}/bin/build_fshd_report.py" \
  "$sample_id" \
  "$classification_dir" \
  "$subset_dir" \
  "$flagstat_txt" \
  "$coverage_tsv" \
  "$hap_tsv" \
  "$methylation_dir" \
  "$methyl_tsv" \
  "$contraction_threshold" \
  "$locus_bed" \
  "$methyl_bed" \
  "${outdir}/${sample_id}.fshd.report.html" \
  "${outdir}/${sample_id}.fshd.report.summary.tsv"

echo "Wrote ${outdir}/${sample_id}.fshd.report.html"
echo "Wrote ${outdir}/${sample_id}.fshd.report.summary.tsv"
