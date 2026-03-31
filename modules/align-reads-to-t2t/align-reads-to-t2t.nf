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

    bam2fq_threads=${params.align_bam2fq_threads}
    map_threads=${params.align_map_threads}
    sort_threads=${params.align_sort_threads}

    # Inputs may be either mapped BAMs or uBAMs with no @SQ targets.
    if samtools view -H "${input_bam}" | grep -q '^@SQ'; then
      samtools quickcheck -v "${input_bam}"
    else
      samtools quickcheck -u -v "${input_bam}"
    fi

    reset_cmd_body="samtools reset --no-PG -x tp,cm,s1,s2,NM,MD,AS,SA,ms,nn,ts,cg,cs,dv,de,rl"
    fastq_cmd_body="samtools fastq -@ \${bam2fq_threads} -n -T MM,ML"

    # Use reset+fastq so re-alignment does not inherit stale mapping tags from pre-aligned BAMs.
    # MM/ML are preserved through FASTQ comments and restored by minimap2 -y.
    \${reset_cmd_body} "${input_bam}" -o - \
      | \${fastq_cmd_body} - \
      | minimap2 -ax ${params.minimap2_preset} --MD -L -Y -y -2 -K ${params.align_minimap_batch} -t \${map_threads} "${t2t_align_ref}" - \
      | samtools sort -@ \${sort_threads} -o "${sample_id}.t2t.bam" -

    samtools index -@ \${sort_threads} "${sample_id}.t2t.bam"
    aligned_records="\$(samtools view -c "${sample_id}.t2t.bam" || true)"
    if [[ "\${aligned_records}" == "0" ]]; then
      echo "No aligned records were produced for ${sample_id} during T2T alignment." >&2
      exit 1
    fi
    """

  stub:
    """
    set -euo pipefail

    touch "${sample_id}.t2t.bam"
    touch "${sample_id}.t2t.bam.bai"
    """
}
