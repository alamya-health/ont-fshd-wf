#!/usr/bin/env nextflow

process ALIGN_READS_TO_HG38 {

  tag "${sample_id}"

  publishDir "${params.output_dir}/align-reads-to-hg38/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(input_bam), path(hg38_align_ref)

  output:
    tuple val(sample_id), path("${sample_id}.hg38.bam"), path("${sample_id}.hg38.bam.bai")

  script:
    """
    set -euo pipefail

    # Inputs may be either mapped BAMs or uBAMs with no @SQ targets.
    if samtools view -H "${input_bam}" | grep -q '^@SQ'; then
      samtools quickcheck -v "${input_bam}"
    else
      samtools quickcheck -u -v "${input_bam}"
    fi

    reset_cmd_body="samtools reset --no-PG -x tp,cm,s1,s2,NM,MD,AS,SA,ms,nn,ts,cg,cs,dv,de,rl"
    fastq_cmd_body="samtools fastq -n -T MM,ML"

    # Use reset+fastq so re-alignment does not inherit stale mapping tags from pre-aligned BAMs.
    # MM/ML are preserved through FASTQ comments and restored by minimap2 -y.
    \${reset_cmd_body} "${input_bam}" -o - \
      | \${fastq_cmd_body} - \
      | minimap2 -ax ${params.minimap2_preset} --MD -L -Y -y -t ${task.cpus} "${hg38_align_ref}" - \
      | samtools sort -@ ${task.cpus} -o "${sample_id}.hg38.bam" -

    samtools index -@ ${task.cpus} "${sample_id}.hg38.bam"
    aligned_records="\$(samtools view -c "${sample_id}.hg38.bam" || true)"
    if [[ "\${aligned_records}" == "0" ]]; then
      echo "No aligned records were produced for ${sample_id} during hg38 alignment." >&2
      exit 1
    fi
    """

  stub:
    """
    set -euo pipefail

    touch "${sample_id}.hg38.bam"
    touch "${sample_id}.hg38.bam.bai"
    """
}
