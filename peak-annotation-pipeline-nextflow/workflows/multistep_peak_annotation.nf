/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WorkflowMultistepPeakAnnotation } from '../lib/WorkflowMultistepPeakAnnotation'
include { NfcoreSchema } from '../lib/NfcoreSchema'
include { NfcoreTemplate } from '../lib/NfcoreTemplate'
include { GenomeSpeciesMapping } from '../lib/GenomeSpeciesMapping'

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMultistepPeakAnnotation.initialise(params, log)

// Check mandatory parameters - FIXED LOGIC
if (params.input) {
    ch_input = file(params.input)
} else if (!params.consensus_peaks) {
    exit 1, 'Either --input samplesheet or --consensus_peaks must be specified!'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { MACS2_CONSENSUS     } from '../subworkflows/local/macs2_consensus'
include { PEAK_ANNOTATION     } from '../subworkflows/local/peak_annotation'

//
// MODULE: Local modules
//
include { DOWNLOAD_GTF        } from '../modules/local/download_gtf/main'
include { EXTRACT_INTRONS     } from '../modules/local/extract_introns/main'
include { EXTRACT_LNCRNA_MIRNA } from '../modules/local/extract_lncrna_mirna/main'
include { EXPAND_TARGETS_LNCRNA } from '../modules/local/expand_targets_lncrna/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { GUNZIP                      } from '../modules/nf-core/gunzip/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MULTISTEP_PEAK_ANNOTATION {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files (only if --input provided)
    //
    if (params.input) {
        INPUT_CHECK (
            ch_input
        )
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    }

    //
    // Handle GTF file - either provided or download from Ensembl
    //
    ch_gtf = Channel.empty()
    ch_gtf_info = Channel.empty()

    if (params.gtf) {
        // Use provided GTF file
        ch_gtf = Channel.fromPath(params.gtf, checkIfExists: true)
        ch_gtf_info = Channel.value("User provided GTF: ${params.gtf}")
        log.info "Using provided GTF file: ${params.gtf}"
    } else if (params.genome && params.auto_download_references) {
        // Determine species from genome
        def species = params.species ?: GenomeSpeciesMapping.getSpeciesFromGenome(params.genome)
        if (!species) {
            exit 1, "ERROR: Could not determine species for genome '${params.genome}'. Please provide --species or --gtf parameter."
        }

        log.info "Downloading GTF for ${species} from Ensembl release ${params.ensembl_version}"

        // Download GTF from Ensembl
        DOWNLOAD_GTF (
            species,
            params.ensembl_version
        )

        // Decompress GTF if needed
        GUNZIP (
            DOWNLOAD_GTF.out.gtf.map { gtf -> [[id: "${species}_${params.ensembl_version}"], gtf] }
        )

        ch_gtf = GUNZIP.out.gunzip.map { meta, gtf -> gtf }
        ch_gtf_info = DOWNLOAD_GTF.out.info.map { info -> "Downloaded from Ensembl: ${info}" }
        ch_versions = ch_versions.mix(DOWNLOAD_GTF.out.versions)
        ch_versions = ch_versions.mix(GUNZIP.out.versions)
    } else {
        exit 1, "ERROR: Either --gtf must be provided or --genome must be specified with --auto_download_references true"
    }

    //
    // Create channels for other reference files
    //
    ch_crm_bed = params.crm_bed ? Channel.fromPath(params.crm_bed, checkIfExists: true) : Channel.empty()

    // Handle intron file - either provided or extract from GTF
    ch_intron_bed = Channel.empty()
    ch_intron_info = Channel.empty()

    if (params.intron_bed) {
        // Use provided intron file
        ch_intron_bed = Channel.fromPath(params.intron_bed, checkIfExists: true)
        ch_intron_info = Channel.value("User provided intron file: ${params.intron_bed}")
        log.info "Using provided intron file: ${params.intron_bed}"
    } else if (!params.skip_intron && params.auto_download_references) {
        // Extract introns from GTF
        log.info "Extracting first introns from GTF file"

        EXTRACT_INTRONS (
            ch_gtf
        )

        ch_intron_bed = EXTRACT_INTRONS.out.introns
        ch_intron_info = EXTRACT_INTRONS.out.log.map { log -> "Extracted from GTF: ${log}" }
        ch_versions = ch_versions.mix(EXTRACT_INTRONS.out.versions)
    } else if (params.skip_intron) {
        log.info "Skipping intron annotation step"
        ch_intron_info = Channel.value("Intron annotation skipped")
    } else {
        log.info "No intron file provided and auto-extraction disabled"
        ch_intron_info = Channel.value("No intron annotation performed")
    }

    //
    // Extract lncRNA-miRNA relationships from GTF
    //
    ch_lncrna_mirna_mapping = Channel.empty()
    ch_lncrna_mirna_info = Channel.empty()

    if (params.enable_lncrna_mirna_expansion) {
        log.info "Extracting lncRNA-miRNA relationships from GTF file"

        EXTRACT_LNCRNA_MIRNA (
            ch_gtf
        )

        ch_lncrna_mirna_mapping = EXTRACT_LNCRNA_MIRNA.out.mapping
        ch_lncrna_mirna_info = EXTRACT_LNCRNA_MIRNA.out.log.map { log -> "Extracted lncRNA-miRNA relationships: ${log}" }
        ch_versions = ch_versions.mix(EXTRACT_LNCRNA_MIRNA.out.versions)
    } else {
        log.info "lncRNA-miRNA relationship expansion disabled"
        ch_lncrna_mirna_info = Channel.value("lncRNA-miRNA expansion disabled")
    }

    //
    // Print sample processing information (only if using samplesheet)
    //
    if (params.input) {
        INPUT_CHECK.out.peaks
            .map { meta, peaks -> meta.id }
            .unique()
            .collect()
            .subscribe { sample_list ->
                log.info "Processing ${sample_list.size()} samples: ${sample_list.join(', ')}"
            }
    }

    //
    // Print GTF, intron, and lncRNA-miRNA information
    //
    ch_gtf_info.subscribe { info ->
        log.info "GTF source: ${info}"
    }

    ch_intron_info.subscribe { info ->
        log.info "Intron source: ${info}"
    }

    ch_lncrna_mirna_info.subscribe { info ->
        log.info "lncRNA-miRNA: ${info}"
    }

    //
    // CONSENSUS PEAK LOGIC - FIXED to handle three scenarios
    //
    ch_consensus_peaks = Channel.empty()
    ch_consensus_log = Channel.empty()

    if (params.consensus_peaks) {
        // Scenario 1: User provided pre-computed consensus peaks
        // Skip consensus calling entirely and use provided peaks
        ch_consensus_peaks = Channel.fromPath(params.consensus_peaks, checkIfExists: true)
            .map { peaks ->
                def meta = [id: file(peaks).getBaseName()]
                return [meta, peaks]
            }
        log.info "Using provided consensus peaks: ${params.consensus_peaks}"

    } else if (params.input && !params.skip_consensus) {
        // Scenario 2: Generate consensus peaks from samplesheet
        // Run MACS2_CONSENSUS on input peaks
        MACS2_CONSENSUS (
            INPUT_CHECK.out.peaks,
            params.min_consensus_reps
        )
        ch_consensus_peaks = MACS2_CONSENSUS.out.consensus_peaks
        ch_consensus_log = MACS2_CONSENSUS.out.log
        ch_versions = ch_versions.mix(MACS2_CONSENSUS.out.versions)

        // Log consensus peak results
        ch_consensus_log.subscribe { meta, log_file ->
            log.info "Sample '${meta.id}': Consensus peaks generated"
        }

    } else if (params.input && params.skip_consensus) {
        // Scenario 3: Skip consensus entirely, use individual peaks
        // Pass through individual peaks from samplesheet
        ch_consensus_peaks = INPUT_CHECK.out.peaks
        log.info "Skipping consensus calling - using individual peaks for annotation"

    } else {
        exit 1, "ERROR: No valid peak input configuration. Use --input with optional --skip_consensus, or --consensus_peaks"
    }

    //
    // SUBWORKFLOW: Multi-step peak annotation - processes each sample independently
    //
    PEAK_ANNOTATION (
        ch_consensus_peaks,
        ch_crm_bed,
        ch_intron_bed,
        ch_gtf,
        params.homer_distance,
        params.intersect_overlap_fraction,
        params.skip_crm,
        params.skip_intron,
        params.genome,
        ch_lncrna_mirna_mapping,
        params.enable_lncrna_mirna_expansion
    )
    ch_versions = ch_versions.mix(PEAK_ANNOTATION.out.versions)

    // Log annotation results
    PEAK_ANNOTATION.out.final_report.subscribe { meta, report_file ->
        log.info "Sample '${meta.id}': Multi-step annotation completed"
    }

    // Log expansion results
    if (params.enable_lncrna_mirna_expansion) {
        PEAK_ANNOTATION.out.expansion_log.subscribe { meta, log_file ->
            log.info "Sample '${meta.id}': lncRNA-miRNA expansion completed"
        }
    }

    //
    // MODULE: Collate software versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC - combines results from all samples
    //
    workflow_summary    = WorkflowMultistepPeakAnnotation.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMultistepPeakAnnotation.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(PEAK_ANNOTATION.out.multiqc_files.collect{it[1]}.ifEmpty([]))

    if (ch_consensus_log && !params.skip_consensus && !params.consensus_peaks) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_consensus_log.collect{it[1]}.ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    //
    // Print final summary
    //
    workflow.onComplete {
        log.info "Pipeline completed successfully!"
        if (params.input) {
            def sample_count = 0
            INPUT_CHECK.out.peaks
                .map { meta, peaks -> meta.id }
                .unique()
                .count()
                .subscribe { count ->
                    sample_count = count
                    log.info "Processed ${sample_count} samples with multi-step peak annotation"
                }
        } else {
            log.info "Processed consensus peaks with multi-step peak annotation"
        }

        if (params.gtf) {
            log.info "GTF file used: ${params.gtf}"
        } else {
            log.info "GTF downloaded from Ensembl release ${params.ensembl_version}"
        }
        log.info "Results available in: ${params.outdir}"
    }

    emit:
    multiqc_report = MULTIQC.out.report
    versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// END OF SCRIPT
