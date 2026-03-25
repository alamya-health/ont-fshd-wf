#!/usr/bin/env nextflow

process CALL_FSHD_VARIANTS {

  tag "${sample_id}"

  publishDir "${params.output_dir}/call-fshd-variants/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(hg38_bam), path(hg38_bai), path(t2t_bam), path(t2t_bai)
    path hg38_ref_fasta
    path clair3_model_asset
    path clinvar_vcf_gz
    path clinvar_vcf_tbi
    path snpeff_data_tgz

  output:
    tuple val(sample_id), path("${sample_id}.variant-calling")
    tuple val(sample_id), path("${sample_id}.variant.summary.tsv")

  script:
    """
    set -euo pipefail

    if [[ ! -f "${hg38_bam}.bai" ]]; then
      ln -sf "${hg38_bai}" "${hg38_bam}.bai"
    fi
    if [[ ! -f "${t2t_bam}.bai" ]]; then
      ln -sf "${t2t_bai}" "${t2t_bam}.bai"
    fi

    outdir="${sample_id}.variant-calling"
    mkdir -p "\${outdir}/clair3" "\${outdir}/refs" "\${outdir}/assets"

    model_dir="/opt/models/${params.clair3_model_name}"
    if [[ ! -f "\${model_dir}/pileup.pt" || ! -f "\${model_dir}/full_alignment.pt" ]]; then
      echo "Bundled Clair3 model not found or incomplete: \${model_dir}" >&2
      exit 1
    fi

    if [[ "\$(basename "${clair3_model_asset}")" != "clair3_model.placeholder" ]]; then
      candidate_model_dir=""
      if [[ -d "${clair3_model_asset}" ]]; then
        candidate_model_dir="\$(find "${clair3_model_asset}" -mindepth 0 -maxdepth 4 -type f -name pileup.pt | sort | head -n 1 | xargs -r dirname)"
      elif [[ "${clair3_model_asset}" == *.tar.gz ]]; then
        tar -xzf "${clair3_model_asset}" -C "\${outdir}/assets"
        candidate_model_dir="\$(find "\${outdir}/assets" -mindepth 1 -maxdepth 4 -type f -name pileup.pt | sort | head -n 1 | xargs -r dirname)"
      else
        echo "Unsupported Clair3 model override asset: ${clair3_model_asset}" >&2
        echo "Provide either a directory containing pileup.pt/full_alignment.pt or a .tar.gz archive of that directory." >&2
        exit 1
      fi

      if [[ -n "\${candidate_model_dir}" && -f "\${candidate_model_dir}/full_alignment.pt" ]]; then
        model_dir="\${candidate_model_dir}"
        echo "Using external Clair3 model override at \${model_dir}" >&2
      else
        echo "External Clair3 model override did not contain pileup.pt and full_alignment.pt." >&2
        echo "Leave the external Clair3 model param blank to use bundled model ${params.clair3_model_name}." >&2
        exit 1
      fi
    else
      echo "Using bundled Clair3 model at \${model_dir}" >&2
    fi

    tar -xzf "${snpeff_data_tgz}" -C "\${outdir}/assets"
    snpeff_data_root="\${outdir}/assets"
    if [[ ! -d "\${snpeff_data_root}/data" ]]; then
      candidate="\$(find "\${outdir}/assets" -mindepth 1 -maxdepth 3 -type d -name data | sort | head -n 1)"
      if [[ -n "\${candidate}" ]]; then
        snpeff_data_root="\$(dirname "\${candidate}")"
      fi
    fi
    if [[ ! -d "\${snpeff_data_root}/data/hg38" ]]; then
      echo "Unable to find snpEff hg38 database under extracted assets" >&2
      exit 1
    fi

    cp "${clinvar_vcf_gz}" "\${outdir}/refs/clinvar.vcf.gz"
    cp "${clinvar_vcf_tbi}" "\${outdir}/refs/clinvar.vcf.gz.tbi"

    total_reads="\$(samtools view -c "${hg38_bam}")"
    if [[ "\${total_reads}" -eq 0 ]]; then
      {
        printf "status\\ttotal_clair3_variants\\trelevant_variant_records\\tclinvar_fshd_hits\\thg38_sv_records\\tt2t_sv_records\\trelevant_genes\\n"
        printf "empty_input\\t0\\t0\\t0\\t0\\t0\\tNA\\n"
      } > "${sample_id}.variant.summary.tsv"
      touch "\${outdir}/empty.txt"
      exit 0
    fi

    run_clair3.sh \
      --bam_fn="${hg38_bam}" \
      --ref_fn="${hg38_ref_fasta}" \
      --threads=${task.cpus} \
      --platform=ont \
      --model_path="\${model_dir}" \
      --output="\${outdir}/clair3" \
      --sample_name="${sample_id}" \
      --enable_phasing \
      --whatshap="\$(command -v whatshap)"

    clair3_vcfgz=""
    if [[ -f "\${outdir}/clair3/phased_merge_output.vcf.gz" ]]; then
      clair3_vcfgz="\${outdir}/clair3/phased_merge_output.vcf.gz"
    elif [[ -f "\${outdir}/clair3/merge_output.vcf.gz" ]]; then
      clair3_vcfgz="\${outdir}/clair3/merge_output.vcf.gz"
    else
      echo "Clair3 did not produce a merged VCF" >&2
      exit 1
    fi

    whatshap haplotag \
      -o "\${outdir}/${sample_id}.hg38.haplotagged.bam" \
      --reference "${hg38_ref_fasta}" \
      "\${clair3_vcfgz}" \
      "${hg38_bam}" \
      --ignore-read-groups \
      --output-haplotag-list "\${outdir}/${sample_id}.haplotag-list.tsv" \
      --output-threads=${task.cpus}

    samtools index -@ ${task.cpus} "\${outdir}/${sample_id}.hg38.haplotagged.bam"

    whatshap stats \
      --gtf="\${outdir}/${sample_id}.haploblocks.gtf" \
      "\${clair3_vcfgz}" \
      > "\${outdir}/${sample_id}.whatshap.stats.txt"

    sniffles \
      --input "\${outdir}/${sample_id}.hg38.haplotagged.bam" \
      --vcf "\${outdir}/${sample_id}.hg38.sniffles2.phased.vcf" \
      --phase \
      --output-rnames

    sniffles \
      --input "${t2t_bam}" \
      --vcf "\${outdir}/${sample_id}.t2t.sniffles2.vcf" \
      --output-rnames

    java -Xmx4g -jar /opt/snpeff/snpEff.jar \
      -dataDir "\${snpeff_data_root}" \
      hg38 \
      -noStats \
      -canon \
      "\${clair3_vcfgz}" \
      > "\${outdir}/${sample_id}.hg38.snpeff.vcf"

    bgzip -f "\${outdir}/${sample_id}.hg38.snpeff.vcf"
    tabix -f -p vcf "\${outdir}/${sample_id}.hg38.snpeff.vcf.gz"

    java -Xmx2g -jar /opt/snpeff/SnpSift.jar \
      annotate -v \
      "\${outdir}/refs/clinvar.vcf.gz" \
      "\${outdir}/${sample_id}.hg38.snpeff.vcf.gz" \
      > "\${outdir}/${sample_id}.hg38.snpeff.clinvar.vcf"

    bgzip -f "\${outdir}/${sample_id}.hg38.snpeff.clinvar.vcf"
    tabix -f -p vcf "\${outdir}/${sample_id}.hg38.snpeff.clinvar.vcf.gz"

    filter_expr="((ANN[*].GENE = 'DUX4') | (ANN[*].GENE = 'SMCHD1') | (ANN[*].GENE = 'LRIF1') | (ANN[*].GENE = 'DNMT3B') | (ANN[*].GENE = 'TRIM43') | (ANN[*].GENE = 'CAPN3') | (ANN[*].GENE = 'VCP') | (CLNDN =~ 'Facioscapulohumeral'))"

    java -Xmx2g -jar /opt/snpeff/SnpSift.jar \
      filter "\${filter_expr}" \
      "\${outdir}/${sample_id}.hg38.snpeff.clinvar.vcf.gz" \
      > "\${outdir}/${sample_id}.fshd_relevant.vcf"

    bgzip -f "\${outdir}/${sample_id}.fshd_relevant.vcf"
    tabix -f -p vcf "\${outdir}/${sample_id}.fshd_relevant.vcf.gz"

    write_fshd_variant_summary.py \
      "${sample_id}" \
      "\${clair3_vcfgz}" \
      "\${outdir}/${sample_id}.fshd_relevant.vcf.gz" \
      "\${outdir}/${sample_id}.hg38.sniffles2.phased.vcf" \
      "\${outdir}/${sample_id}.t2t.sniffles2.vcf" \
      "${sample_id}.variant.summary.tsv"
    """

  stub:
    """
    set -euo pipefail

    mkdir -p "${sample_id}.variant-calling"
    touch "${sample_id}.variant-calling/${sample_id}.hg38.haplotagged.bam"
    touch "${sample_id}.variant-calling/${sample_id}.hg38.haplotagged.bam.bai"
    touch "${sample_id}.variant-calling/${sample_id}.hg38.sniffles2.phased.vcf"
    touch "${sample_id}.variant-calling/${sample_id}.t2t.sniffles2.vcf"
    cat <<'EOF' > "${sample_id}.variant.summary.tsv"
status	total_clair3_variants	relevant_variant_records	clinvar_fshd_hits	hg38_sv_records	t2t_sv_records	relevant_genes
ready	1	1	0	0	0	SMCHD1
EOF
    """
}
