#!/usr/bin/env nextflow

process MERGE_UNALIGNED_BAMS {

  tag "${sample_id}"

  publishDir "${params.output_dir}/merge-unaligned-bams/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(input_ubam_dir)

  output:
    tuple val(sample_id), path("${sample_id}.input.bam")

  script:
    def mergeThreads = Math.max((task.cpus as int) - 1, 1)

    """
    set -euo pipefail

    bam_root="${input_ubam_dir}"
    declare -a bam_files=()

    if [[ -f "\$bam_root" ]]; then
      case "\$bam_root" in
        *.bam|*.BAM)
          bam_files=("\$bam_root")
          ;;
        *)
          echo "Input path is a file but not a BAM: \$bam_root" >&2
          exit 1
          ;;
      esac
    elif [[ -d "\$bam_root" ]]; then
      mapfile -d '' -t bam_files < <(find -L "\$bam_root" -type f \\( -iname '*.bam' \\) -print0 | sort -z)
    else
      echo "Input path does not exist or is not accessible: \$bam_root" >&2
      exit 1
    fi

    if [[ "\${#bam_files[@]}" -eq 0 ]]; then
      echo "No BAM files found in \$bam_root" >&2
      exit 1
    fi

    if [[ "\${#bam_files[@]}" -eq 1 ]]; then
      cp "\${bam_files[0]}" "${sample_id}.input.bam"
    else
      samtools merge -u -@ ${mergeThreads} -f "${sample_id}.input.bam" "\${bam_files[@]}"
    fi

    # uBAMs may not carry @SQ targets; validate them in unmapped mode.
    samtools quickcheck -u -v "${sample_id}.input.bam"
    """

  stub:
    """
    set -euo pipefail

    touch "${sample_id}.input.bam"
    """
}
