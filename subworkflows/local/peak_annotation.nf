// subworkflows/local/peak_annotation.nf

//
// Multi-step peak annotation subworkflow
//

include { BEDTOOLS_INTERSECT_CRM         } from '../../modules/local/bedtools_intersect_crm/main'
include { BEDTOOLS_INTERSECT_INTRON      } from '../../modules/local/bedtools_intersect_intron/main'
include { HOMER_ANNOTATEPEAKS                                  } from '../../modules/nf-core/homer/annotatepeaks/main'
include { PLOT_HOMER_ANNOTATEPEAKS                             } from '../../modules/local/plot_homer_annotatepeaks/main'
include { PREPARE_ANNOTATION_OUTPUT                            } from '../../modules/local/prepare_annotation_output/main'
include { EXPAND_TARGETS_LNCRNA                                } from '../../modules/local/expand_targets_lncrna/main'

workflow PEAK_ANNOTATION {

    take:
    consensus_peaks           // channel: [ val(meta), path(peaks) ]
    crm_bed                   // path: CRM BED file (value channel)
    intron_bed                // path: intron BED file (value channel)
    gtf                       // path: GTF file (value channel)
    homer_distance            // val: HOMER distance parameter
    overlap_fraction          // val: overlap fraction for intersections
    skip_crm                  // val: skip CRM annotation
    skip_intron               // val: skip intron annotation
    genome                    // val: genome name
    lncrna_mirna_mapping      // path: lncRNA-miRNA mapping (optional)
    enable_lncrna_expansion   // val: enable lncRNA-miRNA expansion

    main:

    ch_versions = Channel.empty()

    // Initialize channels for tracking annotation outputs
    ch_crm_annotation = Channel.empty()
    ch_intron_annotation = Channel.empty()
    ch_homer_annotation = Channel.empty()

    // Track sample metadata
    ch_sample_meta = consensus_peaks.map { meta, peaks -> meta }

    //
    // Step 1: CRM annotation (optional)
    //
    ch_peaks_for_intron = consensus_peaks

    if (!skip_crm && crm_bed) {
        log.info "Running CRM annotation step"

        BEDTOOLS_INTERSECT_CRM (
            consensus_peaks,
            crm_bed,
            overlap_fraction
        )
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_CRM.out.versions)

        // Get peaks that don't intersect with CRM for next step
        ch_peaks_for_intron = BEDTOOLS_INTERSECT_CRM.out.non_intersected
        ch_crm_annotation = BEDTOOLS_INTERSECT_CRM.out.annotation

        log.info "CRM annotation completed"
    } else {
        log.info "Skipping CRM annotation step"
        // Create empty CRM annotation files for samples
        ch_crm_annotation = ch_sample_meta.map { meta ->
            [meta, file("NO_FILE_CRM")]
        }
    }

    //
    // Step 2: Intron annotation (optional)
    //
    ch_peaks_for_homer = ch_peaks_for_intron

    if (!skip_intron && intron_bed) {
        log.info "Running intron annotation step"

        BEDTOOLS_INTERSECT_INTRON (
            ch_peaks_for_intron,
            intron_bed,
            overlap_fraction
        )
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_INTRON.out.versions)

        // Get peaks that don't intersect with introns for HOMER
        ch_peaks_for_homer = BEDTOOLS_INTERSECT_INTRON.out.non_intersected
        ch_intron_annotation = BEDTOOLS_INTERSECT_INTRON.out.annotation

        log.info "Intron annotation completed"
    } else {
        log.info "Skipping intron annotation step"
        // Create empty intron annotation files for samples
        ch_intron_annotation = ch_sample_meta.map { meta ->
            [meta, file("NO_FILE_INTRON")]
        }
    }

    //
    // Step 3: HOMER annotation (always run on remaining peaks)
    //
    log.info "Running HOMER annotation step"

    HOMER_ANNOTATEPEAKS (
        ch_peaks_for_homer,
        gtf,
        genome,
        homer_distance
    )
    ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions)
    ch_homer_annotation = HOMER_ANNOTATEPEAKS.out.txt

    //
    // Step 4: Plot HOMER results (optional)
    //
    if (!params.skip_plots) {
        PLOT_HOMER_ANNOTATEPEAKS (
            HOMER_ANNOTATEPEAKS.out.txt
        )
        ch_versions = ch_versions.mix(PLOT_HOMER_ANNOTATEPEAKS.out.versions)
    }

    //
    // Step 5: PREPARE FINAL ANNOTATION OUTPUT - THIS WAS MISSING!
    //
    log.info "Preparing final annotation output"

    // Combine all annotation results by sample - join by metadata
    ch_all_annotations = ch_sample_meta
        .join(ch_crm_annotation, by: 0, remainder: true)
        .join(ch_intron_annotation, by: 0, remainder: true)
        .join(ch_homer_annotation, by: 0, remainder: true)
        .map { meta, crm_file, intron_file, homer_file ->
            [meta, [crm_file, intron_file, homer_file]]
        }

    // Get the original consensus peaks for the final output process
    ch_consensus_for_output = consensus_peaks

    PREPARE_ANNOTATION_OUTPUT (
        ch_all_annotations,
        ch_consensus_for_output
    )
    ch_versions = ch_versions.mix(PREPARE_ANNOTATION_OUTPUT.out.versions)

    //
    // Step 6: lncRNA-miRNA expansion (optional)
    //
    ch_expansion_log = Channel.empty()
    ch_expanded_targets = PREPARE_ANNOTATION_OUTPUT.out.all_genes

    if (enable_lncrna_expansion && lncrna_mirna_mapping) {
        log.info "Running lncRNA-miRNA target expansion"

        EXPAND_TARGETS_LNCRNA (
            PREPARE_ANNOTATION_OUTPUT.out.all_genes,
            lncrna_mirna_mapping
        )
        ch_versions = ch_versions.mix(EXPAND_TARGETS_LNCRNA.out.versions)
        ch_expansion_log = EXPAND_TARGETS_LNCRNA.out.log
        ch_expanded_targets = EXPAND_TARGETS_LNCRNA.out.expanded_targets

        log.info "lncRNA-miRNA expansion completed"
    } else {
        log.info "Skipping lncRNA-miRNA expansion"
    }

    emit:
    target_genes       = ch_expanded_targets                              // channel: [ val(meta), path(target_genes) ]
    final_report       = PREPARE_ANNOTATION_OUTPUT.out.report             // channel: [ val(meta), path(report) ]
    summary_log        = PREPARE_ANNOTATION_OUTPUT.out.summary            // channel: [ val(meta), path(log) ]
    expansion_log      = ch_expansion_log                                 // channel: [ val(meta), path(log) ]

    // Individual annotation outputs for MultiQC
    crm_annotation     = ch_crm_annotation                                // channel: [ val(meta), path(annotation) ]
    intron_annotation  = ch_intron_annotation                             // channel: [ val(meta), path(annotation) ]
    homer_annotation   = ch_homer_annotation                              // channel: [ val(meta), path(annotation) ]
    homer_plots        = params.skip_plots ? Channel.empty() : PLOT_HOMER_ANNOTATEPEAKS.out.plots

    // MultiQC files
    multiqc_files      = Channel.empty()
        .mix(PREPARE_ANNOTATION_OUTPUT.out.multiqc_files.map{it[1]})
        .mix(ch_expansion_log.map{it[1]})
        .mix(params.skip_plots ? Channel.empty() : PLOT_HOMER_ANNOTATEPEAKS.out.plots.map{it[1]})

    versions           = ch_versions                                      // channel: [ path(versions) ]
}
