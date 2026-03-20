#!/usr/bin/env nextflow

process PAS_CHECK_FSHD {

  tag "${sample_id}"

  publishDir "${params.output_dir}/pas-check-fshd/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(locus_bam), path(locus_bai)

  output:
    tuple val(sample_id), path("${sample_id}.PAS.txt")

  script:
    """
    set -euo pipefail

    if [[ ! -f "${locus_bam}.bai" ]]; then
      ln -sf "${locus_bai}" "${locus_bam}.bai"
    fi

    check_fshd_pas.py "${locus_bam}" "${sample_id}.PAS.txt"
    """

  stub:
    """
    cat <<'EOF' > "${sample_id}.PAS.txt"
read_id	pas_sequence	pas_type	pas_status
stub_read	ATTAAAAT	disrupted	stub
EOF
    """
}
