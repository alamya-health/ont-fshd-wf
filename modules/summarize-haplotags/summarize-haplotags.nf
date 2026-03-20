#!/usr/bin/env nextflow

process SUMMARIZE_HAPLOTAGS {

  tag "${sample_id}"

  publishDir "${params.output_dir}/summarize-haplotags/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(original_bam), path(classification_dir)

  output:
    tuple val(sample_id), path("${sample_id}.haplotag.summary.tsv")

  script:
    """
    set -euo pipefail

    summarize_fshd_haplotags.py "${original_bam}" "${classification_dir}" "${sample_id}.haplotag.summary.tsv"
    """

  stub:
    """
    cat <<'EOF' > "${sample_id}.haplotag.summary.tsv"
sample_id	haplotype	read_count
stub	4qA	1
EOF
    """
}
