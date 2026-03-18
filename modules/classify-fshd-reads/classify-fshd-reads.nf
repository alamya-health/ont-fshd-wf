#!/usr/bin/env nextflow

process CLASSIFY_FSHD_READS {

  tag "${sample_id}"

  publishDir "${params.output_dir}/classify-fshd-reads/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(blast_txt), path(pas_txt)
    path fshd_analysis_rscript

  output:
    tuple val(sample_id), path("${sample_id}.fshd.classification")

  script:
    """
    set -euo pipefail

    mkdir "${sample_id}.fshd.classification"

    Rscript "${fshd_analysis_rscript}" "${blast_txt}" "${pas_txt}" .

    shopt -s nullglob
    for f in *.csv *.txt *.log; do
      mv "$f" "${sample_id}.fshd.classification/"
    done
    if [[ -d blast_results ]]; then
      mv blast_results "${sample_id}.fshd.classification/"
    fi
    """
}
