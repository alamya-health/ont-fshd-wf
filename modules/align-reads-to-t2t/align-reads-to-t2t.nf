#!/usr/bin/env nextflow

process ALIGN_READS_TO_T2T {

  tag "${sample_id}"

  publishDir "${params.output_dir}/align-reads-to-t2t/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(input_bam), path(t2t_align_ref)

  output:
    tuple val(sample_id), path("${sample_id}.t2t.bam"), path("${sample_id}.t2t.bam.bai")

  script:
    """
    set -euo pipefail

    # Inputs may be either mapped BAMs or uBAMs with no @SQ targets.
    if samtools view -H "${input_bam}" | grep -q '^@SQ'; then
      samtools quickcheck -v "${input_bam}"
    else
      samtools quickcheck -u -v "${input_bam}"
    fi

    samtools fastq \
      -F 0x900 \
      -n \
      -T MM,ML \
      -0 /dev/null \
      -s /dev/null \
      "${input_bam}" \
      | minimap2 -ax ${params.minimap2_preset} --MD -L -Y -y -t ${task.cpus} "${t2t_align_ref}" - \
      | samtools sort -@ ${task.cpus} -o "${sample_id}.t2t.bam" -

    samtools index -@ ${task.cpus} "${sample_id}.t2t.bam"
    """

  stub:
    """
    set -euo pipefail

    touch "${sample_id}.t2t.bam"
    touch "${sample_id}.t2t.bam.bai"
    """
}
