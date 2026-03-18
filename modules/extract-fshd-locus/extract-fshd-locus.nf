#!/usr/bin/env nextflow

process EXTRACT_FSHD_LOCUS {

  tag "${sample_id}"

  publishDir "${params.output_dir}/extract-fshd-locus/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(input_bam), path(input_bai), path(target_bed)

  output:
    tuple val(sample_id), path("${sample_id}.fshd.locus.bam"), path("${sample_id}.fshd.locus.bam.bai")
    tuple val(sample_id), path("${sample_id}.fshd.locus.flagstat.txt")
    tuple val(sample_id), path("${sample_id}.fshd.locus.coverage.tsv")

  script:
    """
    set -euo pipefail

    if [[ ! -f "${input_bam}.bai" ]]; then
      ln -sf "${input_bai}" "${input_bam}.bai"
    fi

    samtools quickcheck -v "${input_bam}"

    samtools view \
      -@ ${task.cpus} \
      -h \
      -b \
      -L "${target_bed}" \
      -o "${sample_id}.fshd.locus.bam" \
      "${input_bam}"

    samtools index -@ ${task.cpus} "${sample_id}.fshd.locus.bam"

    samtools flagstat -@ ${task.cpus} "${sample_id}.fshd.locus.bam" > "${sample_id}.fshd.locus.flagstat.txt"
    samtools coverage -b "${target_bed}" "${sample_id}.fshd.locus.bam" > "${sample_id}.fshd.locus.coverage.tsv"
    """
}
