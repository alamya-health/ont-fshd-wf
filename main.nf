#!/usr/bin/nextflow
nextflow.enable.dsl=2

include { MERGE_UNALIGNED_BAMS } from "./modules/merge-unaligned-bams/merge-unaligned-bams.nf"
include { BUILD_T2T_MMI } from "./modules/build-t2t-mmi/build-t2t-mmi.nf"
include { ALIGN_READS_TO_T2T } from "./modules/align-reads-to-t2t/align-reads-to-t2t.nf"
include { EXTRACT_FSHD_LOCUS } from "./modules/extract-fshd-locus/extract-fshd-locus.nf"
include { BLAST_FSHD_READS } from "./modules/blast-fshd-reads/blast-fshd-reads.nf"
include { PAS_CHECK_FSHD } from "./modules/pas-check-fshd/pas-check-fshd.nf"
include { CLASSIFY_FSHD_READS } from "./modules/classify-fshd-reads/classify-fshd-reads.nf"
include { EXTRACT_CLASSIFIED_READ_SUBSETS } from "./modules/extract-classified-read-subsets/extract-classified-read-subsets.nf"
include { SUMMARIZE_HAPLOTAGS } from "./modules/summarize-haplotags/summarize-haplotags.nf"
include { PROFILE_FSHD_METHYLATION } from "./modules/profile-fshd-methylation/profile-fshd-methylation.nf"
include { BUILD_FSHD_REPORT } from "./modules/build-fshd-report/build-fshd-report.nf"
include { SUMMARIZE_COHORT } from "./modules/summarize-cohort/summarize-cohort.nf"

workflow {

    main:

    def usingSamplesheet = params.samplesheet ? params.samplesheet.toString().trim() : ''
    def projectBlastDb = file(params.blast_db_dir.toString())
    def projectLocusBed = file(params.fshd_locus_bed.toString())
    def projectMethylBed = file(params.methylation_regions_bed.toString())
    def projectMethylPileupBed = file(params.methylation_pileup_bed.toString())
    def projectRScript = file(params.fshd_analysis_rscript.toString())
    def getString = { row, key, fallback = '' ->
        (row[key] ?: fallback ?: '').toString().trim()
    }
    def getBoolean = { row, key, fallback = false ->
        def v = row[key]
        if( v == null || v.toString().trim() == '' ) {
            v = fallback
        }
        if( v instanceof Boolean ) {
            return v
        }
        return v.toString().trim().toLowerCase() in ['1', 'true', 't', 'yes', 'y']
    }

    ch_samples_meta = null
    if( usingSamplesheet ) {
        ch_samples_meta = Channel
            .fromPath(params.samplesheet.toString())
            .ifEmpty { error "Samplesheet not found or empty: ${params.samplesheet}" }
            .splitCsv(header: true)
            .map { row ->
                def sid = getString(row, 'sample_id')
                if( !sid ) {
                    throw new IllegalArgumentException("Samplesheet row missing required column: sample_id")
                }
                tuple(sid, row)
            }
    }
    else {
        if( !params.sample_id ) {
            error "Missing required parameter: provide either --samplesheet or --sample_id"
        }
        ch_samples_meta = Channel.value(tuple(params.sample_id.toString(), [:]))
    }

    if( !params.report_only ) {
        ch_samples_meta = ch_samples_meta.map { sid, row ->
            def ubamDir = getString(row, 'input_ubam_dir', params.input_ubam_dir)
            def alignedBam = getString(row, 'input_aligned_bam', params.input_aligned_bam)
            if( ubamDir && alignedBam ) {
                throw new IllegalArgumentException("Provide either input_ubam_dir or input_aligned_bam for sample '${sid}', not both")
            }
            if( !ubamDir && !alignedBam ) {
                throw new IllegalArgumentException("Missing input for sample '${sid}': provide input_ubam_dir or input_aligned_bam")
            }
            tuple(sid, row)
        }
    }

    ch_ubam_samples = Channel.empty()
    ch_aligned_samples = Channel.empty()
    ch_aligned_reuse = Channel.empty()
    ch_aligned_realign = Channel.empty()

    if( !params.report_only ) {
        ch_ubam_samples = ch_samples_meta
            .filter { sid, row ->
                def ubamDir = getString(row, 'input_ubam_dir', params.input_ubam_dir)
                return ubamDir ? true : false
            }
            .map { sid, row ->
                tuple(sid, file(getString(row, 'input_ubam_dir', params.input_ubam_dir)))
            }

        ch_aligned_samples = ch_samples_meta
            .filter { sid, row ->
                def ubamDir = getString(row, 'input_ubam_dir', params.input_ubam_dir)
                def alignedBam = getString(row, 'input_aligned_bam', params.input_aligned_bam)
                return !ubamDir && alignedBam
            }

        ch_aligned_reuse = ch_aligned_samples
            .filter { sid, row ->
                getBoolean(row, 'reuse_input_t2t_alignment', params.reuse_input_t2t_alignment ?: false)
            }
            .map { sid, row ->
                def bamPath = getString(row, 'input_aligned_bam', params.input_aligned_bam)
                def baiPath = getString(row, 'input_aligned_bai', params.input_aligned_bai)
                if( !baiPath ) {
                    throw new IllegalArgumentException("Aligned input reuse requires input_aligned_bai for sample '${sid}'")
                }
                tuple(sid, file(bamPath), file(baiPath))
            }

        ch_aligned_realign = ch_aligned_samples
            .filter { sid, row ->
                !getBoolean(row, 'reuse_input_t2t_alignment', params.reuse_input_t2t_alignment ?: false)
            }
            .map { sid, row ->
                tuple(sid, file(getString(row, 'input_aligned_bam', params.input_aligned_bam)))
            }
    }

    if( !params.output_dir ) {
        error "Missing required parameter: --output_dir"
    }

    if( params.report_only ) {
        if( !params.report_only_input_root ) {
            error "Missing required parameter for report-only mode: --report_only_input_root"
        }
    }
    else {
        if( !params.t2t_ref_fasta ) {
            error "Missing required parameter: --t2t_ref_fasta"
        }
    }

    def ch_report_html = null
    def ch_report_summary = null

    if( params.report_only ) {
        def reportRoot = params.report_only_input_root.toString().trim().replaceAll('/+$', '')
        ch_report_static = ch_samples_meta.map { sid, row ->
            def classification_dir = file("${reportRoot}/classify-fshd-reads/${sid}/${sid}.fshd.classification")
            def subset_dir = file("${reportRoot}/extract-classified-read-subsets/${sid}/${sid}.classified-subsets")
            def flagstat_txt = file("${reportRoot}/extract-fshd-locus/${sid}/${sid}.fshd.locus.flagstat.txt")
            def coverage_tsv = file("${reportRoot}/extract-fshd-locus/${sid}/${sid}.fshd.locus.coverage.tsv")
            def haplotag_summary_tsv = file("${reportRoot}/summarize-haplotags/${sid}/${sid}.haplotag.summary.tsv")
            def methylation_dir = file("${reportRoot}/profile-fshd-methylation/${sid}/${sid}.methylation")
            def methylation_summary_tsv = file("${reportRoot}/profile-fshd-methylation/${sid}/${sid}.methylation.summary.tsv")
            tuple(sid, classification_dir, subset_dir, flagstat_txt, coverage_tsv, haplotag_summary_tsv, methylation_dir, methylation_summary_tsv)
        }

        if( params.report_only_run_methylation ) {
            ch_report_methylation_input = ch_samples_meta.map { sid, row ->
                def subset_dir = file("${reportRoot}/extract-classified-read-subsets/${sid}/${sid}.classified-subsets")
                def explicitOriginal = getString(row, 'report_only_original_bam', params.report_only_original_bam)
                def originalBam = explicitOriginal ? file(explicitOriginal) : file("${reportRoot}/merge-unaligned-bams/${sid}/${sid}.input.bam")
                if( !originalBam.exists() ) {
                    throw new IllegalArgumentException("Report-only methylation rerun requires an original donor BAM for sample '${sid}'. Set --report_only_original_bam or ensure ${reportRoot}/merge-unaligned-bams/${sid}/${sid}.input.bam exists.")
                }
                tuple(sid, originalBam, subset_dir)
            }

            PROFILE_FSHD_METHYLATION(
                ch_report_methylation_input,
                Channel.value(file(params.t2t_ref_fasta.toString())),
                Channel.value(projectMethylBed),
                Channel.value(projectMethylPileupBed)
            )

            ch_report_input = ch_report_static
                .map { sid, classification_dir, subset_dir, flagstat_txt, coverage_tsv, haplotag_summary_tsv, methylation_dir, methylation_summary_tsv ->
                    tuple(sid, classification_dir, subset_dir, flagstat_txt, coverage_tsv, haplotag_summary_tsv)
                }
                .join(PROFILE_FSHD_METHYLATION.out[0])
                .join(PROFILE_FSHD_METHYLATION.out[1])
                .map { sid, classification_dir, subset_dir, flagstat_txt, coverage_tsv, haplotag_summary_tsv, methylation_dir, methylation_summary_tsv ->
                    tuple(sid, classification_dir, subset_dir, flagstat_txt, coverage_tsv, haplotag_summary_tsv, methylation_dir, methylation_summary_tsv)
                }
        }
        else {
            ch_report_input = ch_report_static
        }

        BUILD_FSHD_REPORT( ch_report_input, Channel.value(projectLocusBed), Channel.value(projectMethylBed) )

        ch_report_html = BUILD_FSHD_REPORT.out[0]
        ch_report_summary = BUILD_FSHD_REPORT.out[1]
    }
    else {

        MERGE_UNALIGNED_BAMS( ch_ubam_samples )
        ch_merged_original_bam = MERGE_UNALIGNED_BAMS.out[0]

        ch_t2t_align_ref = null
        if( params.t2t_ref_mmi ) {
            ch_t2t_align_ref = Channel.value(file(params.t2t_ref_mmi.toString()))
        }
        else {
            BUILD_T2T_MMI(Channel.value(file(params.t2t_ref_fasta.toString())))
            ch_t2t_align_ref = BUILD_T2T_MMI.out
        }

        ch_need_alignment = ch_merged_original_bam.mix(ch_aligned_realign)
        ch_align_input = ch_need_alignment
            .combine(ch_t2t_align_ref)
            .map { sid, bam, t2tAlignRef ->
                tuple(sid, bam, t2tAlignRef)
            }

        ALIGN_READS_TO_T2T( ch_align_input )
        ch_aligned_t2t_bam = ALIGN_READS_TO_T2T.out[0]

        ch_original_bam = ch_merged_original_bam
            .mix(ch_aligned_realign)
            .mix(ch_aligned_reuse.map { sid, bam, bai -> tuple(sid, bam) })

        ch_analysis_bam = ch_aligned_t2t_bam
            .mix(ch_aligned_reuse)

        ch_extract_input = ch_analysis_bam.map { sid, bam, bai ->
            tuple(sid, bam, bai, projectLocusBed)
        }
        EXTRACT_FSHD_LOCUS( ch_extract_input )

        ch_locus_bam = EXTRACT_FSHD_LOCUS.out[0]
        ch_locus_flagstat = EXTRACT_FSHD_LOCUS.out[1]
        ch_locus_coverage = EXTRACT_FSHD_LOCUS.out[2]

        ch_blast_input = ch_locus_bam.map { sid, bam, bai ->
            tuple(sid, bam, projectBlastDb)
        }
        BLAST_FSHD_READS( ch_blast_input )

        PAS_CHECK_FSHD( ch_locus_bam )

        ch_classify_input = BLAST_FSHD_READS.out[0]
            .join(PAS_CHECK_FSHD.out[0])
            .map { sid, blast_txt, pas_txt ->
                tuple(sid, blast_txt, pas_txt)
            }

        CLASSIFY_FSHD_READS(
            ch_classify_input,
            Channel.value(projectRScript)
        )
        ch_classification_dir = CLASSIFY_FSHD_READS.out[0]

        ch_subset_input = ch_locus_bam
            .join(ch_classification_dir)
            .map { sid, locus_bam, locus_bai, classification_dir ->
                tuple(sid, locus_bam, locus_bai, classification_dir)
            }

        EXTRACT_CLASSIFIED_READ_SUBSETS( ch_subset_input )
        ch_subset_dir = EXTRACT_CLASSIFIED_READ_SUBSETS.out[0]

        ch_haplotag_input = ch_original_bam
            .join(ch_classification_dir)
            .map { sid, original_bam, classification_dir ->
                tuple(sid, original_bam, classification_dir)
            }

        SUMMARIZE_HAPLOTAGS( ch_haplotag_input )
        ch_haplotag_summary = SUMMARIZE_HAPLOTAGS.out[0]

        ch_methylation_input = ch_original_bam
            .join(ch_subset_dir)
            .map { sid, original_bam, subset_dir ->
                tuple(sid, original_bam, subset_dir)
            }

        PROFILE_FSHD_METHYLATION(
            ch_methylation_input,
            Channel.value(file(params.t2t_ref_fasta.toString())),
            Channel.value(projectMethylBed),
            Channel.value(projectMethylPileupBed)
        )
        ch_methylation_dir = PROFILE_FSHD_METHYLATION.out[0]
        ch_methylation_summary = PROFILE_FSHD_METHYLATION.out[1]

        ch_report_input = ch_classification_dir
            .join(ch_subset_dir)
            .join(ch_locus_flagstat)
            .join(ch_locus_coverage)
            .join(ch_haplotag_summary)
            .join(ch_methylation_dir)
            .join(ch_methylation_summary)
            .map { sid, classification_dir, subset_dir, flagstat_txt, coverage_tsv, haplotag_summary_tsv, methylation_dir, methylation_summary_tsv ->
                tuple(sid, classification_dir, subset_dir, flagstat_txt, coverage_tsv, haplotag_summary_tsv, methylation_dir, methylation_summary_tsv)
            }

        BUILD_FSHD_REPORT( ch_report_input, Channel.value(projectLocusBed), Channel.value(projectMethylBed) )

        ch_report_html = BUILD_FSHD_REPORT.out[0]
        ch_report_summary = BUILD_FSHD_REPORT.out[1]
    }

    ch_report_summary_files = ch_report_summary.map { sid, summary_tsv -> summary_tsv }.collect()
    SUMMARIZE_COHORT( ch_report_summary_files )

    emit:
    report_html = ch_report_html
    report_summary = ch_report_summary
    cohort_summary_tsv = SUMMARIZE_COHORT.out[0]
    cohort_summary_html = SUMMARIZE_COHORT.out[1]
}
