//
// Multi-step peak annotation subworkflow
//

include { BEDTOOLS_INTERSECT_CRM     } from '../../modules/local/bedtools_intersect_crm/main'
include { BEDTOOLS_INTERSECT_INTRON  } from '../../modules/local/bedtools_intersect_intron/main'
include { HOMER_ANNOTATEPEAKS        } from '../../modules/nf-core/homer/annotatepeaks/main'
include { PREPARE_ANNOTATION_OUTPUT  } from '../../modules/local/prepare_annotation_output/main'
include { EXPAND_TARGETS_LNCRNA      } from '../../modules/local/expand_targets_lncrna/main'

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
            overlap_fraction
        )
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_CRM.out.versions)
        ch_peaks_for_intron = BEDTOOLS_INTERSECT_CRM.out.non_intersected
        ch_crm_annotation = BEDTOOLS_INTERSECT_CRM.out.annotation
    } else {
        ch_crm_annotation = ch_sample_meta.map { meta -> [meta, file("NO_FILE_CRM")] }
    }

    // Step 2: Intron annotation
    ch_peaks_for_homer = ch_peaks_for_intron

    if (!skip_intron && intron_bed) {
        BEDTOOLS_INTERSECT_INTRON (
            ch_peaks_for_intron,
            intron_bed,
            overlap_fraction
        )
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT_INTRON.out.versions)
        ch_peaks_for_homer = BEDTOOLS_INTERSECT_INTRON.out.non_intersected
        ch_intron_annotation = BEDTOOLS_INTERSECT_INTRON.out.annotation
    } else {
        ch_intron_annotation = ch_sample_meta.map { meta -> [meta, file("NO_FILE_INTRON")] }
    }

    // Step 3: HOMER annotation
    HOMER_ANNOTATEPEAKS (
        ch_peaks_for_homer,   // tuple val(meta), path(peak)
        file(genome),         // path fasta
        file(gtf)            // path gtf
    )
    ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions)
    ch_homer_annotation = HOMER_ANNOTATEPEAKS.out.txt

    // Step 4: Final output
    ch_all_annotations = ch_sample_meta
        .join(ch_crm_annotation, by: 0, remainder: true)
        .join(ch_intron_annotation, by: 0, remainder: true)
        .join(ch_homer_annotation, by: 0, remainder: true)
        .map { meta, crm_file, intron_file, homer_file ->
            [meta, [crm_file, intron_file, homer_file]]
        }

    PREPARE_ANNOTATION_OUTPUT (
        ch_all_annotations,
        consensus_peaks
    )
    ch_versions = ch_versions.mix(PREPARE_ANNOTATION_OUTPUT.out.versions)

    // Step 5: lncRNA expansion
    ch_expanded_targets = PREPARE_ANNOTATION_OUTPUT.out.all_genes
    ch_expansion_log = Channel.empty()

    if (enable_lncrna_expansion && lncrna_mirna_mapping) {
        EXPAND_TARGETS_LNCRNA (
            PREPARE_ANNOTATION_OUTPUT.out.all_genes,
            lncrna_mirna_mapping
        )
        ch_versions = ch_versions.mix(EXPAND_TARGETS_LNCRNA.out.versions)
        ch_expansion_log = EXPAND_TARGETS_LNCRNA.out.log
        ch_expanded_targets = EXPAND_TARGETS_LNCRNA.out.expanded_targets
    }

    emit:
    target_genes       = ch_expanded_targets
    final_report       = PREPARE_ANNOTATION_OUTPUT.out.report
    summary_log        = PREPARE_ANNOTATION_OUTPUT.out.summary
    expansion_log      = ch_expansion_log
    crm_annotation     = ch_crm_annotation
    intron_annotation  = ch_intron_annotation
    homer_annotation   = ch_homer_annotation
    homer_plots        = Channel.empty()
    multiqc_files      = PREPARE_ANNOTATION_OUTPUT.out.multiqc_files.map{it[1]}
    versions           = ch_versions
}
