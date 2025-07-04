//
// This file holds several functions specific to the workflow/multistep_peak_annotation.nf in the nf-core/multisteppeak-annotation pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowMultistepPeakAnnotation {

    //
    // Check and validate parameters - FIXED VERSION
    //
    public static void initialise(params, log) {

        // Check input requirements
        if (!params.input && !params.consensus_peaks) {
            Nextflow.error("Either --input (samplesheet) or --consensus_peaks must be specified!")
        }

        // Check mutually exclusive parameters
        if (params.input && params.consensus_peaks) {
            Nextflow.error("Cannot specify both --input and --consensus_peaks. Use --input for consensus calling or --consensus_peaks for pre-computed peaks.")
        }

        // If using consensus_peaks, skip_consensus should be implied
        if (params.consensus_peaks && !params.skip_consensus) {
            log.info "Using provided consensus peaks - consensus calling will be skipped"
        }

        // Check genome parameter
        if (!params.genome && !params.gtf && !params.auto_download_references) {
            Nextflow.error("Either --genome or --gtf must be specified, or --auto_download_references must be enabled!")
        }

        // Check consensus parameters (only relevant if generating consensus)
        if (params.input && !params.skip_consensus) {
            if (params.min_consensus_reps < 1) {
                Nextflow.error("--min_consensus_reps must be >= 1")
            }
        }

        // Check CRM file if annotation is enabled
        if (!params.skip_crm && !params.crm_bed) {
            log.warn("CRM annotation is enabled but no CRM BED file provided via --crm_bed. CRM annotation will be skipped.")
        }

        // Check HOMER distance parameter
        if (params.homer_distance < 0) {
            Nextflow.error("--homer_distance must be >= 0")
        }

        // Check overlap fraction parameter
        if (params.intersect_overlap_fraction < 0 || params.intersect_overlap_fraction > 1) {
            Nextflow.error("--intersect_overlap_fraction must be between 0 and 1")
        }

        // Check Ensembl version parameter
        if (params.ensembl_version < 1) {
            Nextflow.error("--ensembl_version must be >= 1")
        }

        // Validate lncRNA-miRNA parameters
        if (params.enable_lncrna_mirna_expansion && !params.auto_download_references && !params.gtf) {
            Nextflow.error("lncRNA-miRNA expansion requires either --gtf or --auto_download_references to be enabled")
        }

        // Print parameter summary
        paramsSummaryLog(params, log)
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"

        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //
    public static String methodsDescriptionText(run_workflow, mqc_methods_description_text) {
        // Convert to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta.manifest = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta.doi_text = meta.manifest.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest.doi}\'>${meta.manifest.doi}</a>)" : ""
        meta.nodoi_text = meta.manifest.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta.tool_citations = toolCitationText().replaceAll(", \\.", ".")
        meta.tool_bibliography = toolBibliographyText()

        def methods_text = mqc_methods_description_text.text

        def engine = new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Generate tool citations for MultiQC
    //
    public static String toolCitationText() {
        def citation_text = [
            "Tools used in the workflow included:",
            "BEDTools (Quinlan and Hall 2010),",
            "HOMER (Heinz et al. 2010),",
            "MACS2 (Zhang et al. 2008),",
            "MultiQC (Ewels et al. 2016)."
        ].join(' ').trim()

        return citation_text
    }

    //
    // Generate tool bibliography for MultiQC
    //
    public static String toolBibliographyText() {
        def bibliography_text = [
            "<li>Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841-842.</li>",
            "<li>Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432</li>",
            "<li>Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS. Model-based analysis of ChIP-Seq (MACS). Genome Biol. 2008;9(9):R137. doi: 10.1186/gb-2008-9-9-r137. Epub 2008 Sep 17. PMID: 18798982; PMCID: PMC2592715.</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

        return bibliography_text
    }

    //
    // Print parameter summary log to screen
    //
    public static void paramsSummaryLog(params, log) {
        Map colors = NfcoreTemplate.logColours(params.monochrome_logs)
        log.info NfcoreTemplate.dashedLine(params.monochrome_logs)
        log.info "${colors.green}Multi-step Peak Annotation Parameters${colors.reset}"
        log.info "${colors.dim}  input                           : ${colors.reset} ${params.input}"
        log.info "${colors.dim}  consensus_peaks                 : ${colors.reset} ${params.consensus_peaks ?: 'None'}"
        log.info "${colors.dim}  genome                          : ${colors.reset} ${params.genome ?: 'None'}"
        log.info "${colors.dim}  gtf                             : ${colors.reset} ${params.gtf ?: 'None'}"
        log.info "${colors.dim}  crm_bed                         : ${colors.reset} ${params.crm_bed ?: 'None'}"
        log.info "${colors.dim}  intron_bed                      : ${colors.reset} ${params.intron_bed ?: 'None'}"
        log.info "${colors.dim}  min_consensus_reps              : ${colors.reset} ${params.min_consensus_reps}"
        log.info "${colors.dim}  homer_distance                  : ${colors.reset} ${params.homer_distance}"
        log.info "${colors.dim}  intersect_overlap_fraction      : ${colors.reset} ${params.intersect_overlap_fraction}"
        log.info "${colors.dim}  auto_download_references        : ${colors.reset} ${params.auto_download_references}"
        log.info "${colors.dim}  ensembl_version                 : ${colors.reset} ${params.ensembl_version}"
        log.info "${colors.dim}  species                         : ${colors.reset} ${params.species ?: 'Auto-detect'}"
        log.info "${colors.dim}  enable_lncrna_mirna_expansion   : ${colors.reset} ${params.enable_lncrna_mirna_expansion}"
        log.info "${colors.dim}  skip_consensus                  : ${colors.reset} ${params.skip_consensus}"
        log.info "${colors.dim}  skip_crm                        : ${colors.reset} ${params.skip_crm}"
        log.info "${colors.dim}  skip_intron                     : ${colors.reset} ${params.skip_intron}"
        log.info NfcoreTemplate.dashedLine(params.monochrome_logs)
    }

    //
    // Validate peak calling parameters
    //
    public static void validatePeakParams(params, log) {
        // Check that input files exist
        if (params.input && !file(params.input).exists()) {
            Nextflow.error("Input samplesheet does not exist: ${params.input}")
        }

        if (params.gtf && !file(params.gtf).exists()) {
            Nextflow.error("GTF file does not exist: ${params.gtf}")
        }

        if (params.crm_bed && !file(params.crm_bed).exists()) {
            Nextflow.error("CRM BED file does not exist: ${params.crm_bed}")
        }

        if (params.intron_bed && !file(params.intron_bed).exists()) {
            Nextflow.error("Intron BED file does not exist: ${params.intron_bed}")
        }

        if (params.consensus_peaks && !file(params.consensus_peaks).exists()) {
            Nextflow.error("Consensus peaks file does not exist: ${params.consensus_peaks}")
        }

        // Validate output directory
        if (!params.outdir) {
            Nextflow.error("Output directory not specified!")
        }

        log.info "Parameter validation completed successfully"
    }

    //
    // Get attribute from genome config file e.g. fasta
    //
    public static Object getGenomeAttribute(params, attribute) {
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[params.genome].containsKey(attribute)) {
                return params.genomes[params.genome][attribute]
            }
        }
        return null
    }

    //
    // Create summary of pipeline parameters
    //
    public static LinkedHashMap paramsSummaryMap(workflow, params) {
        // Core Nextflow options
        def summary = [:]
        summary['Pipeline Name']         = workflow.manifest.name
        summary['Pipeline Version']      = workflow.manifest.version
        summary['Run Name']              = workflow.runName
        summary['User']                  = workflow.userName
        summary['Profile']               = workflow.profile
        summary['Config Files']          = workflow.configFiles
        summary['Container']             = workflow.container
        summary['Container Engine']      = workflow.containerEngine
        summary['Launch Dir']            = workflow.launchDir
        summary['Working Dir']           = workflow.workDir
        summary['Project Dir']           = workflow.projectDir
        summary['Script Dir']            = workflow.projectDir
        summary['Command Line']          = workflow.commandLine
        summary['Nextflow Version']      = workflow.nextflow.version
        summary['Nextflow Build']        = workflow.nextflow.build
        summary['Nextflow Compile Time'] = workflow.nextflow.timestamp

        // Input/output options
        summary['Input']                 = params.input
        summary['Consensus Peaks']       = params.consensus_peaks ?: 'None'
        summary['Output dir']            = params.outdir
        summary['Tracedir']              = params.tracedir
        summary['Publish dir mode']      = params.publish_dir_mode
        summary['MultiQC config']        = params.multiqc_config ?: 'None'
        summary['MultiQC title']         = params.multiqc_title ?: 'None'
        summary['Max Memory']            = params.max_memory
        summary['Max CPUs']              = params.max_cpus
        summary['Max Time']              = params.max_time

        // Pipeline specific options
        summary['Genome']                = params.genome ?: 'None'
        summary['GTF']                   = params.gtf ?: 'None'
        summary['CRM BED']               = params.crm_bed ?: 'None'
        summary['Intron BED']            = params.intron_bed ?: 'None'
        summary['Min Consensus Reps']    = params.min_consensus_reps
        summary['HOMER Distance']        = params.homer_distance
        summary['Intersect Overlap']     = params.intersect_overlap_fraction
        summary['Auto Download Refs']    = params.auto_download_references
        summary['Ensembl Version']       = params.ensembl_version
        summary['Species']               = params.species ?: 'Auto-detect'
        summary['lncRNA-miRNA Expansion'] = params.enable_lncrna_mirna_expansion
        summary['Skip Consensus']        = params.skip_consensus
        summary['Skip CRM']              = params.skip_crm
        summary['Skip Intron']           = params.skip_intron

        return summary
    }
}

// END OF SCRIPT
