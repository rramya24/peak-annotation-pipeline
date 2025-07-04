//
// Peak annotation subworkflow
//

include { BEDTOOLS_INTERSECT_CRM    } from '../../modules/local/bedtools_intersect_crm/main'
include { BEDTOOLS_INTERSECT_INTRON } from '../../modules/local/bedtools_intersect_intron/main'
include { HOMER_ANNOTATEPEAKS       } from '../../modules/nf-core/homer/annotatepeaks/main'
include { CONVERT_GENEID_TO_SYMBOL  } from '../../modules/local/convert_geneid_to_symbol/main'
include { PREPARE_ANNOTATION_OUTPUT } from '../../modules/local/prepare_annotation_output/main'
include { EXPAND_TARGETS_LNCRNA     } from '../../modules/local/expand_targets_lncrna/main'

workflow PEAK_ANNOTATION {
    take:
    consensus_peaks         // channel: [ val(meta), path(consensus_peaks) ]
    crm_bed                 // channel: path(crm_bed)
    intron_bed              // channel: path(intron_bed)
    gtf                     // channel: path(gtf)
    homer_distance          // integer: distance for HOMER annotation
    intersect_overlap       // float: overlap fraction for bedtools intersect
    skip_crm                // boolean: skip CRM annotation
    skip_intron            // boolean: skip intron annotation
    genome                  // string: genome name for HOMER
    lncrna_mirna_mapping    // channel: path(lncrna_mirna_mapping)
    enable_lncrna_expansion // boolean: enable lncRNA-miRNA expansion

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Initialize channels for different annotation steps
    ch_crm_annotated = Channel.empty()
    ch_intron_annotated = Channel.empty()
    ch_homer_annotated = Channel.empty()
    ch_remaining_peaks = consensus_peaks

    //
    // STEP 1: CRM annotation (if not skipped)
    //
    if (!skip_crm && crm_bed) {
        BEDTOOLS_INTERSECT_CRM (
            ch_remaining_peaks,
            crm_bed,
            intersect_overlap
        )
        ch_crm_annotated = BEDTOOLS_INTERSECT_CRM.out.intersected
        ch_remaining_peaks = BEDTOOLS_INTERSECT_CRM.out.non_intersected
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_CRM.out.versions)
    }

    //
    // STEP 2: Intron annotation (if not skipped)
    //
    if (!skip_intron && intron_bed) {
        BEDTOOLS_INTERSECT_INTRON (
            ch_remaining_peaks,
            intron_bed,
            intersect_overlap
        )
        ch_intron_annotated = BEDTOOLS_INTERSECT_INTRON.out.intersected
        ch_remaining_peaks = BEDTOOLS_INTERSECT_INTRON.out.non_intersected
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_INTRON.out.versions)
    }

    //
    // STEP 3: HOMER annotation for remaining peaks
    //
    if (ch_remaining_peaks) {
        HOMER_ANNOTATEPEAKS (
            ch_remaining_peaks,
            gtf,
            genome,
            homer_distance
        )
        ch_homer_annotated = HOMER_ANNOTATEPEAKS.out.txt
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions)
    }

    //
    // STEP 4: Convert Gene IDs to Gene Symbols
    //
    ch_all_annotations = Channel.empty()
    ch_all_annotations = ch_all_annotations.mix(ch_crm_annotated.map{ meta, file -> [meta, file, 'crm'] })
    ch_all_annotations = ch_all_annotations.mix(ch_intron_annotated.map{ meta, file -> [meta, file, 'intron'] })
    ch_all_annotations = ch_all_annotations.mix(ch_homer_annotated.map{ meta, file -> [meta, file, 'homer'] })

    CONVERT_GENEID_TO_SYMBOL (
        ch_all_annotations,
        gtf
    )
    ch_versions = ch_versions.mix(CONVERT_GENEID_TO_SYMBOL.out.versions)

    //
    // STEP 5: Prepare final annotation output
    //
    PREPARE_ANNOTATION_OUTPUT (
        CONVERT_GENEID_TO_SYMBOL.out.gene_symbols.groupTuple(),
        consensus_peaks
    )
    ch_versions = ch_versions.mix(PREPARE_ANNOTATION_OUTPUT.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_ANNOTATION_OUTPUT.out.multiqc_files)

    //
    // STEP 6: Expand targets with lncRNA-miRNA relationships (if enabled)
    //
    ch_final_targets = PREPARE_ANNOTATION_OUTPUT.out.all_genes
    ch_expansion_log = Channel.empty()

    if (enable_lncrna_expansion && lncrna_mirna_mapping) {
        EXPAND_TARGETS_LNCRNA (
            PREPARE_ANNOTATION_OUTPUT.out.all_genes,
            lncrna_mirna_mapping
        )
        ch_final_targets = EXPAND_TARGETS_LNCRNA.out.expanded_targets
        ch_expansion_log = EXPAND_TARGETS_LNCRNA.out.log
        ch_versions = ch_versions.mix(EXPAND_TARGETS_LNCRNA.out.versions)
    }

    emit:
    crm_annotated     = ch_crm_annotated                              // channel: [ val(meta), path(crm_annotated) ]
    intron_annotated  = ch_intron_annotated                           // channel: [ val(meta), path(intron_annotated) ]
    homer_annotated   = ch_homer_annotated                            // channel: [ val(meta), path(homer_annotated) ]
    gene_symbols      = CONVERT_GENEID_TO_SYMBOL.out.gene_symbols     // channel: [ val(meta), path(gene_symbols) ]
    final_report      = PREPARE_ANNOTATION_OUTPUT.out.report          // channel: [ val(meta), path(report) ]
    final_targets     = ch_final_targets                              // channel: [ val(meta), path(final_targets) ]
    expansion_log     = ch_expansion_log                              // channel: [ val(meta), path(expansion_log) ]
    multiqc_files     = ch_multiqc_files                              // channel: [ val(meta), path(multiqc_files) ]
    versions          = ch_versions                                   // channel: [ versions.yml ]
}

# END OF SUBWORKFLOW
