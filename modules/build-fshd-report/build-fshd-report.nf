#!/usr/bin/env nextflow

process BUILD_FSHD_REPORT {

  tag "${sample_id}"

  publishDir "${params.output_dir}/build-fshd-report/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(classification_dir), path(subset_dir), path(flagstat_txt), path(coverage_tsv), path(haplotag_summary_tsv), path(methylation_summary_tsv), path(variant_summary_tsv)

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
      "${methylation_summary_tsv}" \
      "${variant_summary_tsv}" \
      "${params.contraction_threshold_ru}" \
      "${sample_id}.fshd.report.html" \
      "${sample_id}.fshd.report.summary.tsv"
    """
}
