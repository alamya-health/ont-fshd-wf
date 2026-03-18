#!/usr/bin/env nextflow

process PROFILE_FSHD_METHYLATION {

  tag "${sample_id}"

  publishDir "${params.output_dir}/profile-fshd-methylation/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(original_bam), path(subset_dir)
    path t2t_ref_fasta
    path methylation_bed

  output:
    tuple val(sample_id), path("${sample_id}.methylation")
    tuple val(sample_id), path("${sample_id}.methylation.summary.tsv")

  script:
    """
    set -euo pipefail

    outdir="${sample_id}.methylation"
    mkdir -p "\${outdir}"
    summary="${sample_id}.methylation.summary.tsv"
    printf "subset_name\\tstatus\\tsubset_read_count\\tdonor_read_count\\tmean_percent_modified\\tpileup_bed_gz\\tstats_tsv\\trepaired_bam\\n" > "\${summary}"

    if [[ "${params.run_methylation}" != "true" ]]; then
      for subset_name in 4qA_all 4qA_complete chimeric; do
        printf "%s\\t%s\\t0\\t0\\tNA\\tNA\\tNA\\tNA\\n" "\${subset_name}" "disabled" >> "\${summary}"
      done
      exit 0
    fi

    if python3 - "${original_bam}" <<'PY'
import sys
import pysam

bam_path = sys.argv[1]
bam = pysam.AlignmentFile(bam_path, "rb", check_sq=False)
found = False
for idx, aln in enumerate(bam.fetch(until_eof=True)):
    if aln.has_tag("MM") or aln.has_tag("Mm"):
        found = True
        break
    if idx >= 2000:
        break
bam.close()
sys.exit(0 if found else 1)
PY
    then
      modtags_ok=0
    else
      modtags_ok=1
    fi

    if [[ "\${modtags_ok}" -ne 0 ]]; then
      for subset_name in 4qA_all 4qA_complete chimeric; do
        printf "%s\\t%s\\t0\\t0\\tNA\\tNA\\tNA\\tNA\\n" "\${subset_name}" "no_mod_tags_detected" >> "\${summary}"
      done
      exit 0
    fi

    declare -a subsets=(
      "4qA_all|${sample_id}.4qA_all.bam"
      "4qA_complete|${sample_id}.4qA_complete.bam"
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

      samtools view "\${subset_bam}" | cut -f1 | sort -u > "\${work_prefix}.ids.txt"
      samtools view -F 0x900 -N "\${work_prefix}.ids.txt" -b "${original_bam}" > "\${work_prefix}.donor.raw.bam"
      donor_read_count="\$(samtools view -c "\${work_prefix}.donor.raw.bam" || true)"
      if [[ "\${donor_read_count}" == "0" ]]; then
        printf "%s\\t%s\\t%s\\t0\\tNA\\tNA\\tNA\\tNA\\n" "\${subset_name}" "no_matching_donor_reads" "\${subset_read_count}" >> "\${summary}"
        continue
      fi

      samtools sort -n -o "\${work_prefix}.donor.name.bam" "\${work_prefix}.donor.raw.bam"
      samtools sort -n -o "\${work_prefix}.acceptor.name.bam" "\${subset_bam}"

      modkit repair \
        "\${work_prefix}.donor.name.bam" \
        "\${work_prefix}.acceptor.name.bam" \
        "\${work_prefix}.repaired.name.bam"

      samtools sort -o "\${work_prefix}.repaired.bam" "\${work_prefix}.repaired.name.bam"
      samtools index "\${work_prefix}.repaired.bam"

      modkit pileup \
        "\${work_prefix}.repaired.bam" \
        "\${work_prefix}.methyl.bed" \
        --cpg \
        --ref "${t2t_ref_fasta}"

      bgzip -f "\${work_prefix}.methyl.bed"
      tabix -f "\${work_prefix}.methyl.bed.gz"

      modkit stats \
        --regions "${methylation_bed}" \
        -o "\${work_prefix}.stats.tsv" \
        "\${work_prefix}.methyl.bed.gz"

      mean_pct="\$(python3 - "\${work_prefix}.methyl.bed.gz" <<'PY'
import gzip
import sys

bed_path = sys.argv[1]
weighted = 0.0
coverage = 0.0
with gzip.open(bed_path, "rt") as fh:
    for line in fh:
        if not line.strip() or line.startswith("track"):
            continue
        parts = line.rstrip("\\n").split("\\t")
        if len(parts) < 11:
            continue
        try:
            cov = float(parts[9])
            pct = float(parts[10])
        except ValueError:
            continue
        weighted += cov * pct
        coverage += cov
print(f"{(weighted / coverage) if coverage else 0.0:.3f}")
PY
)"

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
}
