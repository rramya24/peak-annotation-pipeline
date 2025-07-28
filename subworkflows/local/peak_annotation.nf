//
// Multi-step peak annotation subworkflow
//

include { BEDTOOLS_INTERSECT_CRM     } from '../../modules/local/bedtools_intersect_crm/main'
include { BEDTOOLS_INTERSECT_INTRON  } from '../../modules/local/bedtools_intersect_intron/main'
include { HOMER_ANNOTATEPEAKS        } from '../../modules/nf-core/homer/annotatepeaks/main'
include { PREPARE_ANNOTATION_OUTPUT  } from '../../modules/local/prepare_annotation_output/main'

workflow PEAK_ANNOTATION {

    take:
    consensus_peaks
    crm_bed
    intron_bed
    gtf
    homer_distance
    overlap_fraction
    skip_crm
    skip_intron
    genome
    lncrna_mirna_mapping
    enable_lncrna_expansion

    main:

    ch_versions = Channel.empty()
    ch_sample_meta = consensus_peaks.map { meta, peaks -> meta }

    // Step 1: CRM annotation
    ch_peaks_for_intron = consensus_peaks

    if (!skip_crm && crm_bed) {
        BEDTOOLS_INTERSECT_CRM (
            consensus_peaks,
            crm_bed,
            gtf,               // Pass GTF for gene name lookup
            overlap_fraction
        )
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_CRM.out.versions)
        ch_peaks_for_intron = BEDTOOLS_INTERSECT_CRM.out.non_intersected
        ch_crm_annotation = BEDTOOLS_INTERSECT_CRM.out.annotation
    } else {
        ch_crm_annotation = ch_sample_meta.map { meta -> [meta, file("NO_FILE_CRM")] }
    }

    // Step 2: Intron annotation (gets non-intersected peaks from CRM step)
    ch_peaks_for_homer = ch_peaks_for_intron

    if (!skip_intron && intron_bed) {
        BEDTOOLS_INTERSECT_INTRON (
            ch_peaks_for_intron,
            intron_bed,
            gtf,               // Pass GTF for gene name lookup
            overlap_fraction
        )
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_INTRON.out.versions)
        ch_peaks_for_homer = BEDTOOLS_INTERSECT_INTRON.out.non_intersected
        ch_intron_annotation = BEDTOOLS_INTERSECT_INTRON.out.annotation
    } else {
        ch_intron_annotation = ch_sample_meta.map { meta -> [meta, file("NO_FILE_INTRON")] }
    }

    // Step 3: HOMER annotation (gets non-intersected peaks from intron step)
    HOMER_ANNOTATEPEAKS (
        ch_peaks_for_homer,   // tuple val(meta), path(peak)
        genome,               // path fasta
        gtf                   // path gtf
    )
    ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions)
    ch_homer_annotation = HOMER_ANNOTATEPEAKS.out.txt

    // Step 4: Prepare all annotation files for final output
    ch_all_annotations = ch_sample_meta
        .join(ch_crm_annotation, by: 0, remainder: true)
        .join(ch_intron_annotation, by: 0, remainder: true)
        .join(ch_homer_annotation, by: 0, remainder: true)
        .map { meta, crm_file, intron_file, homer_file ->
            [meta, [crm_file, intron_file, homer_file]]
        }

    // Step 5: Final annotation output with multiple filtering levels and lncRNA-miRNA expansion
    PREPARE_ANNOTATION_OUTPUT (
        ch_all_annotations,    // tuple val(meta), path(gene_symbol_files)
        consensus_peaks,       // tuple val(meta2), path(consensus_peaks)
        gtf,                   // path(gtf_file) - ADDED
        lncrna_mirna_mapping   // path(lncrna_mirna_mapping), optional: true - ADDED
    )
    ch_versions = ch_versions.mix(PREPARE_ANNOTATION_OUTPUT.out.versions)

    emit:
    // Main target outputs
    final_targets      = PREPARE_ANNOTATION_OUTPUT.out.final_targets
    final_detailed     = PREPARE_ANNOTATION_OUTPUT.out.final_detailed

    // Alternative filtering levels
    raw_targets        = PREPARE_ANNOTATION_OUTPUT.out.raw_targets
    raw_detailed       = PREPARE_ANNOTATION_OUTPUT.out.raw_detailed
    exon_filtered_targets = PREPARE_ANNOTATION_OUTPUT.out.exon_filtered_targets
    exon_filtered_detailed = PREPARE_ANNOTATION_OUTPUT.out.exon_filtered_detailed
    biotype_filtered_targets = PREPARE_ANNOTATION_OUTPUT.out.biotype_filtered_targets
    biotype_filtered_detailed = PREPARE_ANNOTATION_OUTPUT.out.biotype_filtered_detailed

    // Legacy and summary outputs
    target_genes       = PREPARE_ANNOTATION_OUTPUT.out.all_genes  // Legacy compatibility
    final_report       = PREPARE_ANNOTATION_OUTPUT.out.final_targets  // Use main output as report
    summary_log        = PREPARE_ANNOTATION_OUTPUT.out.summary
    expansion_log      = PREPARE_ANNOTATION_OUTPUT.out.expansion_log

    // Individual annotation outputs
    crm_annotation     = ch_crm_annotation
    intron_annotation  = ch_intron_annotation
    homer_annotation   = ch_homer_annotation

    // Other outputs
    homer_plots        = Channel.empty()
    multiqc_files      = PREPARE_ANNOTATION_OUTPUT.out.multiqc_files.map{it[1]}  // FIXED THIS LINE
    versions           = ch_versions
}
