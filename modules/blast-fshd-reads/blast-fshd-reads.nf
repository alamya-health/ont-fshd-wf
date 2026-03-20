#!/usr/bin/env nextflow

process BLAST_FSHD_READS {

  tag "${sample_id}"

  publishDir "${params.output_dir}/blast-fshd-reads/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(locus_bam), path(blast_db_dir)

  output:
    tuple val(sample_id), path("${sample_id}.fshd-blast.txt")

  script:
    """
    set -euo pipefail

    samtools fasta -F 0x900 "${locus_bam}" > "${sample_id}.fshd.locus.fasta"

    blastn \
      -db "${blast_db_dir}/FSHD-blast" \
      -query "${sample_id}.fshd.locus.fasta" \
      -num_threads ${task.cpus} \
      -out "${sample_id}.fshd-blast.txt" \
      -outfmt '6 qseqid sseqid pident slen length mismatch gapopen qstart qend sstart send evalue bitscore'
    """

  stub:
    """
    cat <<'EOF' > "${sample_id}.fshd-blast.txt"
stub_read	chr4_D4Z4-4q35	99.0	1000	1000	0	0	1	1000	1	1000	0.0	2000
EOF
    """
}
