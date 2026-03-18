#!/usr/bin/env nextflow

process ALIGN_READS_TO_T2T {

  tag "${sample_id}"

  publishDir "${params.output_dir}/align-reads-to-t2t/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(input_bam), path(t2t_ref_fasta)

  output:
    tuple val(sample_id), path("${sample_id}.t2t.bam"), path("${sample_id}.t2t.bam.bai")

  script:
    """
    set -euo pipefail

    samtools quickcheck -v "${input_bam}"

    samtools fastq -F 0x900 -T '*' "${input_bam}" \
      | minimap2 -ax lr:hq --MD -L -Y -t ${task.cpus} "${t2t_ref_fasta}" - \
      | samtools sort -@ ${task.cpus} -o "${sample_id}.t2t.bam" -

    samtools index -@ ${task.cpus} "${sample_id}.t2t.bam"
    """
}
