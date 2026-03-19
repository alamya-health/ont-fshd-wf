#!/usr/bin/env nextflow

process ALIGN_READS_TO_HG38 {

  tag "${sample_id}"

  publishDir "${params.output_dir}/align-reads-to-hg38/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(input_bam), path(hg38_ref_fasta)

  output:
    tuple val(sample_id), path("${sample_id}.hg38.bam"), path("${sample_id}.hg38.bam.bai")

  script:
    """
    set -euo pipefail

    samtools quickcheck -v "${input_bam}"

    samtools fastq \
      -F 0x900 \
      -n \
      -T MM,ML \
      -0 /dev/null \
      -s /dev/null \
      "${input_bam}" \
      | minimap2 -ax lr:hq --MD -L -Y -y -t ${task.cpus} "${hg38_ref_fasta}" - \
      | samtools sort -@ ${task.cpus} -o "${sample_id}.hg38.bam" -

    samtools index -@ ${task.cpus} "${sample_id}.hg38.bam"
    """
}
