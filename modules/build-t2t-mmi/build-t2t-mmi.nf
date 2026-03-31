#!/usr/bin/env nextflow

process BUILD_T2T_MMI {

  tag "${t2t_ref_fasta.simpleName}"

  publishDir "${params.output_dir}/reference-index/t2t", mode: 'copy', overwrite: true

  input:
    path t2t_ref_fasta

  output:
    path("${t2t_ref_fasta.simpleName}.mmi")

  script:
    """
    set -euo pipefail

    minimap2 -t ${task.cpus} -d "${t2t_ref_fasta.simpleName}.mmi" "${t2t_ref_fasta}"
    """

  stub:
    """
    set -euo pipefail

    touch "${t2t_ref_fasta.simpleName}.mmi"
    """
}
