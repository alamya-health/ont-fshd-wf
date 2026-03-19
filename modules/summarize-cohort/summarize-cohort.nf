#!/usr/bin/env nextflow

process SUMMARIZE_COHORT {

  publishDir "${params.output_dir}/summarize-cohort", mode: 'copy', overwrite: true

  input:
    path summary_files

  output:
    path("cohort.fshd.summary.tsv")
    path("cohort.fshd.summary.html")

  script:
    """
    set -euo pipefail

    summarize_fshd_cohort.py "cohort.fshd.summary.tsv" "cohort.fshd.summary.html" ${summary_files}
    """
}
