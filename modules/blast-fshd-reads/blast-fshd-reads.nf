#!/usr/bin/env nextflow

process BLAST_FSHD_READS {

  tag "${sample_id}"

  publishDir "${params.output_dir}/blast-fshd-reads/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(locus_bam), path(blast_db_dir)

  output:
    tuple val(sample_id), path("${sample_id}.fshd-blast.txt")
    tuple val(sample_id), path("${sample_id}.fshd-blast.metrics.tsv")
    tuple val(sample_id), path("${sample_id}.fshd.locus.fasta")

  script:
    """
    set -euo pipefail

    locus_primary_reads="\$(samtools view -c -F 0x900 "${locus_bam}" || true)"
    samtools fasta -F 0x900 "${locus_bam}" > "${sample_id}.fshd.locus.fasta"
    locus_fasta_reads="\$(grep -c '^>' "${sample_id}.fshd.locus.fasta" || true)"

    blastn \
      -db "${blast_db_dir}/FSHD-blast" \
      -query "${sample_id}.fshd.locus.fasta" \
      -num_threads ${task.cpus} \
      -out "${sample_id}.fshd-blast.txt" \
      -outfmt '6 qseqid sseqid pident slen length mismatch gapopen qstart qend sstart send evalue bitscore'

    blast_hit_rows="\$(wc -l < "${sample_id}.fshd-blast.txt" | tr -d ' ')"
    {
      printf "metric\tvalue\n"
      printf "locus_primary_reads\t%s\n" "\${locus_primary_reads}"
      printf "locus_fasta_reads\t%s\n" "\${locus_fasta_reads}"
      printf "blast_hit_rows\t%s\n" "\${blast_hit_rows}"
    } > "${sample_id}.fshd-blast.metrics.tsv"
    """

  stub:
    """
    cat <<'EOF' > "${sample_id}.fshd.locus.fasta"
>stub_read
ACGT
EOF
    cat <<'EOF' > "${sample_id}.fshd-blast.txt"
stub_read	chr4_D4Z4-4q35	99.0	1000	1000	0	0	1	1000	1	1000	0.0	2000
EOF
    cat <<'EOF' > "${sample_id}.fshd-blast.metrics.tsv"
metric	value
locus_primary_reads	1
locus_fasta_reads	1
blast_hit_rows	1
EOF
    """
}
