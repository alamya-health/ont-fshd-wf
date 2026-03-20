#!/usr/bin/env nextflow

process CLASSIFY_FSHD_READS {

  tag "${sample_id}"

  publishDir "${params.output_dir}/classify-fshd-reads/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(blast_txt), path(pas_txt)
    path fshd_analysis_rscript

  output:
    tuple val(sample_id), path("${sample_id}.fshd.classification")

  script:
    """
    set -euo pipefail

    mkdir "${sample_id}.fshd.classification"

    Rscript "${fshd_analysis_rscript}" "${blast_txt}" "${pas_txt}" .

    shopt -s nullglob
    for f in *.csv *.txt *.log; do
      mv "$f" "${sample_id}.fshd.classification/"
    done
    if [[ -d blast_results ]]; then
      mv blast_results "${sample_id}.fshd.classification/"
    fi
    """

  stub:
    """
    set -euo pipefail

    outdir="${sample_id}.fshd.classification"
    mkdir -p "\${outdir}"
    for id_file in \
      4qA_all-reads-ID.txt \
      4qA_complete-reads-ID.txt \
      4qB_all-reads-ID.txt \
      4qB_complete-reads-ID.txt \
      chimeric-reads-ID.txt \
      chr10_all-reads-ID.txt \
      chr10_complete-reads-ID.txt \
      chr4-undefined_all-reads-ID.txt \
      D4Z4-only_chr4-reads-ID.txt \
      D4Z4-only_chr10-reads-ID.txt
    do
      printf "stub_read\n" > "\${outdir}/\${id_file}"
    done
    cat <<'EOF' > "\${outdir}/${sample_id}.overview.txt"
stub classification
EOF
    """
}
