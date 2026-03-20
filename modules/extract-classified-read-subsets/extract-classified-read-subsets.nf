#!/usr/bin/env nextflow

process EXTRACT_CLASSIFIED_READ_SUBSETS {

  tag "${sample_id}"

  publishDir "${params.output_dir}/extract-classified-read-subsets/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(locus_bam), path(locus_bai), path(classification_dir)

  output:
    tuple val(sample_id), path("${sample_id}.classified-subsets")

  script:
    """
    set -euo pipefail

    if [[ ! -f "${locus_bam}.bai" ]]; then
      ln -sf "${locus_bai}" "${locus_bam}.bai"
    fi

    outdir="${sample_id}.classified-subsets"
    mkdir -p "\${outdir}"

    manifest="\${outdir}/${sample_id}.classified-subsets.manifest.tsv"
    printf "subset_name\\tid_file\\tbam_path\\tread_count\\tstatus\\n" > "\${manifest}"

    declare -a subsets=(
      "4qA_all|4qA_all-reads-ID.txt"
      "4qA_complete|4qA_complete-reads-ID.txt"
      "4qB_all|4qB_all-reads-ID.txt"
      "4qB_complete|4qB_complete-reads-ID.txt"
      "chimeric|chimeric-reads-ID.txt"
      "chr10_all|chr10_all-reads-ID.txt"
      "chr10_complete|chr10_complete-reads-ID.txt"
      "chr4_undefined|chr4-undefined_all-reads-ID.txt"
      "D4Z4_chr4_only|D4Z4-only_chr4-reads-ID.txt"
      "D4Z4_chr10_only|D4Z4-only_chr10-reads-ID.txt"
    )

    for entry in "\${subsets[@]}"; do
      subset_name="\${entry%%|*}"
      id_file="\${entry#*|}"
      id_path="${classification_dir}/\${id_file}"
      out_bam="\${outdir}/${sample_id}.\${subset_name}.bam"

      if [[ ! -s "\${id_path}" ]]; then
        printf "%s\\t%s\\tNA\\t0\\tmissing_or_empty\\n" "\${subset_name}" "\${id_file}" >> "\${manifest}"
        continue
      fi

      samtools view \
        -@ ${task.cpus} \
        -F 0x904 \
        -N "\${id_path}" \
        -b "${locus_bam}" \
        | samtools sort -@ ${task.cpus} -o "\${out_bam}" -

      samtools index -@ ${task.cpus} "\${out_bam}"
      read_count="\$(samtools view -c "\${out_bam}" || true)"
      printf "%s\\t%s\\t%s\\t%s\\tready\\n" "\${subset_name}" "\${id_file}" "\${out_bam}" "\${read_count}" >> "\${manifest}"
    done
    """

  stub:
    """
    set -euo pipefail

    outdir="${sample_id}.classified-subsets"
    mkdir -p "\${outdir}"
    manifest="\${outdir}/${sample_id}.classified-subsets.manifest.tsv"
    cat <<'EOF' > "\${manifest}"
subset_name	id_file	bam_path	read_count	status
4qA_all	4qA_all-reads-ID.txt	stub	1	ready
4qA_complete	4qA_complete-reads-ID.txt	stub	1	ready
chimeric	chimeric-reads-ID.txt	stub	1	ready
EOF
    touch "\${outdir}/${sample_id}.4qA_all.bam"
    touch "\${outdir}/${sample_id}.4qA_complete.bam"
    touch "\${outdir}/${sample_id}.chimeric.bam"
    """
}
