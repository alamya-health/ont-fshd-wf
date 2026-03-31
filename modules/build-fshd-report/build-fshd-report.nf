#!/usr/bin/env nextflow

process BUILD_FSHD_REPORT {

  tag "${sample_id}"

  publishDir "${params.output_dir}/build-fshd-report/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(classification_dir), path(subset_dir), path(flagstat_txt), path(coverage_tsv), path(haplotag_summary_tsv), path(methylation_dir), path(methylation_summary_tsv)
    path locus_bed
    path methylation_bed

  output:
    tuple val(sample_id), path("${sample_id}.fshd.report.html")
    tuple val(sample_id), path("${sample_id}.fshd.report.summary.tsv")

  script:
    """
    set -euo pipefail

    build_fshd_report.py \
      "${sample_id}" \
      "${classification_dir}" \
      "${subset_dir}" \
      "${flagstat_txt}" \
      "${coverage_tsv}" \
      "${haplotag_summary_tsv}" \
      "${methylation_dir}" \
      "${methylation_summary_tsv}" \
      "${params.contraction_threshold_ru}" \
      "${locus_bed}" \
      "${methylation_bed}" \
      "${sample_id}.fshd.report.html" \
      "${sample_id}.fshd.report.summary.tsv"
    """

  stub:
    """
    cat <<EOF > "${sample_id}.fshd.report.html"
<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8"><title>${sample_id} FSHD Report</title></head>
<body><h1>${sample_id} FSHD Report</h1><p>Stub report output.</p></body>
</html>
EOF

    cat <<'EOF' > "${sample_id}.fshd.report.summary.tsv"
sample_id	contraction_status	haplotype_status	methylation_status
stub	stub	stub	stub
EOF
    """
}
