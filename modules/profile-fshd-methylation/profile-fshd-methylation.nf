#!/usr/bin/env nextflow

process PROFILE_FSHD_METHYLATION {

  tag "${sample_id}"

  publishDir "${params.output_dir}/profile-fshd-methylation/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(original_bam), path(subset_dir)
    path t2t_ref_fasta
    path methylation_summary_bed
    path methylation_pileup_bed

  output:
    tuple val(sample_id), path("${sample_id}.methylation")
    tuple val(sample_id), path("${sample_id}.methylation.summary.tsv")

  script:
    """
    set -euo pipefail

    sort_threads=${params.methylation_sort_threads}
    pileup_threads=${params.methylation_pileup_threads}

    outdir="${sample_id}.methylation"
    mkdir -p "\${outdir}"
    summary="${sample_id}.methylation.summary.tsv"
    printf "subset_name\\tstatus\\tsubset_read_count\\tdonor_read_count\\tmean_percent_modified\\tpileup_bed_gz\\tstats_tsv\\trepaired_bam\\n" > "\${summary}"

    if [[ "${params.run_methylation}" != "true" ]]; then
      for subset_name in 4qA_all 4qA_complete 4qB_all chr10_all D4Z4_chr4_only D4Z4_chr10_only chimeric; do
        printf "%s\\t%s\\t0\\t0\\tNA\\tNA\\tNA\\tNA\\n" "\${subset_name}" "disabled" >> "\${summary}"
      done
      exit 0
    fi

    bam_has_modtags() {
      check_bam_modtags.py "\$1"
    }

    if bam_has_modtags "${original_bam}"; then
      original_modtags_ok=0
    else
      original_modtags_ok=1
    fi

    declare -a subsets=(
      "4qA_all|${sample_id}.4qA_all.bam"
      "4qA_complete|${sample_id}.4qA_complete.bam"
      "4qB_all|${sample_id}.4qB_all.bam"
      "chr10_all|${sample_id}.chr10_all.bam"
      "D4Z4_chr4_only|${sample_id}.D4Z4_chr4_only.bam"
      "D4Z4_chr10_only|${sample_id}.D4Z4_chr10_only.bam"
      "chimeric|${sample_id}.chimeric.bam"
    )

    for entry in "\${subsets[@]}"; do
      subset_name="\${entry%%|*}"
      subset_bam="${subset_dir}/\${entry#*|}"
      work_prefix="\${outdir}/\${subset_name}"

      if [[ ! -s "\${subset_bam}" ]]; then
        printf "%s\\t%s\\t0\\t0\\tNA\\tNA\\tNA\\tNA\\n" "\${subset_name}" "missing_subset_bam" >> "\${summary}"
        continue
      fi

      subset_read_count="\$(samtools view -c "\${subset_bam}" || true)"
      if [[ "\${subset_read_count}" == "0" ]]; then
        printf "%s\\t%s\\t0\\t0\\tNA\\tNA\\tNA\\tNA\\n" "\${subset_name}" "empty_subset_bam" >> "\${summary}"
        continue
      fi

      if bam_has_modtags "\${subset_bam}"; then
        donor_read_count="\${subset_read_count}"
        samtools sort -@ \${sort_threads} -o "\${work_prefix}.repaired.bam" "\${subset_bam}"
      elif [[ "\${original_modtags_ok}" -eq 0 ]]; then
        samtools view "\${subset_bam}" | cut -f1 | sort -u > "\${work_prefix}.ids.txt"
        samtools view -@ \${sort_threads} -F 0x900 -N "\${work_prefix}.ids.txt" -b "${original_bam}" > "\${work_prefix}.donor.raw.bam"
        donor_read_count="\$(samtools view -c "\${work_prefix}.donor.raw.bam" || true)"
        if [[ "\${donor_read_count}" == "0" ]]; then
          printf "%s\\t%s\\t%s\\t0\\tNA\\tNA\\tNA\\tNA\\n" "\${subset_name}" "no_matching_donor_reads" "\${subset_read_count}" >> "\${summary}"
          continue
        fi

        samtools sort -@ \${sort_threads} -n -o "\${work_prefix}.donor.name.bam" "\${work_prefix}.donor.raw.bam"
        samtools sort -@ \${sort_threads} -n -o "\${work_prefix}.acceptor.name.bam" "\${subset_bam}"

        modkit repair \
          "\${work_prefix}.donor.name.bam" \
          "\${work_prefix}.acceptor.name.bam" \
          "\${work_prefix}.repaired.name.bam"

        samtools sort -@ \${sort_threads} -o "\${work_prefix}.repaired.bam" "\${work_prefix}.repaired.name.bam"
      else
        printf "%s\\t%s\\t%s\\t0\\tNA\\tNA\\tNA\\tNA\\n" "\${subset_name}" "no_mod_tags_detected" "\${subset_read_count}" >> "\${summary}"
        continue
      fi

      samtools index -@ \${sort_threads} "\${work_prefix}.repaired.bam"

      modkit pileup \
        "\${work_prefix}.repaired.bam" \
        "\${work_prefix}.methyl.bed.gz" \
        --modified-bases 5mC \
        --combine-strands \
        --cpg \
        --threads \${pileup_threads} \
        --interval-size ${params.methylation_interval_size} \
        --include-bed "${methylation_pileup_bed}" \
        --bgzf \
        --ref "${t2t_ref_fasta}"

      tabix -f -p bed "\${work_prefix}.methyl.bed.gz"

      modkit stats \
        --regions "${methylation_summary_bed}" \
        -o "\${work_prefix}.stats.tsv" \
        "\${work_prefix}.methyl.bed.gz"

      mean_pct="\$(compute_modkit_mean_pct.py "\${work_prefix}.methyl.bed.gz")"

      printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \
        "\${subset_name}" \
        "ready" \
        "\${subset_read_count}" \
        "\${donor_read_count}" \
        "\${mean_pct}" \
        "\${work_prefix}.methyl.bed.gz" \
        "\${work_prefix}.stats.tsv" \
        "\${work_prefix}.repaired.bam" >> "\${summary}"
    done
    """

  stub:
    """
    set -euo pipefail

    outdir="${sample_id}.methylation"
    mkdir -p "\${outdir}"
    touch "\${outdir}/4qA_all.repaired.bam"
    touch "\${outdir}/4qA_complete.repaired.bam"
    touch "\${outdir}/chimeric.repaired.bam"
    cat <<'EOF' > "${sample_id}.methylation.summary.tsv"
subset_name	status	subset_read_count	donor_read_count	mean_percent_modified	pileup_bed_gz	stats_tsv	repaired_bam
4qA_all	ready	1	1	42.0	stub	stub	stub
4qA_complete	ready	1	1	40.0	stub	stub	stub
4qB_all	ready	1	1	39.0	stub	stub	stub
chr10_all	ready	1	1	38.0	stub	stub	stub
D4Z4_chr4_only	ready	1	1	41.0	stub	stub	stub
D4Z4_chr10_only	ready	1	1	37.0	stub	stub	stub
chimeric	ready	1	1	35.0	stub	stub	stub
EOF
    """
}
