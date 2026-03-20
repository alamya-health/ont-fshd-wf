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

  stub:
    """
    cat <<'EOF' > cohort.fshd.summary.tsv
sample_id	contraction_status	haplotype_status	methylation_status	variant_status
stub	stub	stub	stub	stub
EOF
    cat <<'EOF' > cohort.fshd.summary.html
<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>FSHD Cohort Summary</title></head><body><h1>FSHD Cohort Summary</h1></body></html>
EOF
    """
}
