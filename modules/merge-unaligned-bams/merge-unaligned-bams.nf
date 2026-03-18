#!/usr/bin/env nextflow

process MERGE_UNALIGNED_BAMS {

  tag "${sample_id}"

  publishDir "${params.output_dir}/merge-unaligned-bams/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(input_ubam_dir)

  output:
    tuple val(sample_id), path("${sample_id}.input.bam")

  script:
    """
    set -euo pipefail

    mapfile -t bam_files < <(find "${input_ubam_dir}" -maxdepth 1 -type f -name '*.bam' | sort)
    if [[ "\${#bam_files[@]}" -eq 0 ]]; then
      echo "No BAM files found in ${input_ubam_dir}" >&2
      exit 1
    fi

    if [[ "\${#bam_files[@]}" -eq 1 ]]; then
      cp "\${bam_files[0]}" "${sample_id}.input.bam"
    else
      samtools merge -f "${sample_id}.input.bam" "\${bam_files[@]}"
    fi

    samtools quickcheck -v "${sample_id}.input.bam"
    """
}
