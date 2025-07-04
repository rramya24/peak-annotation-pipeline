//
// This file holds several functions specific to the main.nf workflow in the nf-core/multisteppeak-annotation pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The pipeline\n" +
            "  https://doi.org/10.5281/zenodo.XXXXXXX\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print workflow version and exit on --version
        if (params.version) {
            String version_string = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${version_string}"
            System.exit(0)
        }

        // Print citation and exit on --help
        if (params.help) {
            def command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh38 -profile docker"
            log.info NfcoreTemplate.logo(workflow, params.monochrome_logs)
            log.info citation(workflow)
            log.info dashedLine(params.monochrome_logs)
            log.info "USAGE:"
            log.info "  The typical command for running the pipeline is as follows:"
            log.info "  ${command}"
            log.info dashedLine(params.monochrome_logs)
            System.exit(0)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (params.enable_conda) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check input has been provided
        if (!params.input) {
            Nextflow.error("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
        }
    }

    //
    // Get attribute from genome config file e.g. fasta
    //
    public static Object getGenomeAttribute(params, attribute) {
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                return params.genomes[ params.genome ][ attribute ]
            }
        }
        return null
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    public static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }

    //
    // Generate methods description for MultiQC
    //
    public static String toolCitationText(params) {
        def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016),",
            "HOMER (Heinz et al. 2010),",
            "BEDTools (Quinlan and Hall 2010),",
            "MACS2 (Zhang et al. 2008)."
        ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {
        def bibliography_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047-3048. doi: /10.1093/bioinformatics/btw354</li>",
            "<li>Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Mol Cell 2010 May 28;38(4):576-589. PMID: 20513432</li>",
            "<li>Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841-842.</li>",
            "<li>Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS. Model-based analysis of ChIP-Seq (MACS). Genome Biol. 2008;9(9):R137. doi: 10.1186/gb-2008-9-9-r137. Epub 2008 Sep 17. PMID: 18798982; PMCID: PMC2592715.</li>"
        ].join(' ').trim()

        return bibliography_text
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
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta.manifest = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta.doi_text = meta.manifest.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest.doi}\'>${meta.manifest.doi}</a>)" : ""
        meta.nodoi_text = meta.manifest.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta.tool_citations = toolCitationText(meta.workflow).replaceAll(", \\.", ".")
        meta.tool_bibliography = toolBibliographyText(meta.workflow)

        def methods_text = mqc_methods_description_text.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if --genome parameter not provided when required
    //
    public static void genomeRequiredError(params, log) {
        if (!params.genome) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  --genome parameter is required but not provided.\n" +
                "  Please provide a genome reference using --genome parameter.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }

    //
    // Print parameter summary log to screen
    //
    public static void paramsSummaryLog(workflow, params, log) {
        Map colors = NfcoreTemplate.logColours(params.monochrome_logs)
        log.info NfcoreTemplate.logo(workflow, params.monochrome_logs)
        log.info NfcoreTemplate.dashedLine(params.monochrome_logs)
        log.info "${colors.green}Core Nextflow options${colors.reset}"
        if (workflow.revision) {
            log.info "${colors.dim}  revision         : ${colors.reset} ${workflow.revision}"
        }
        log.info "${colors.dim}  runName          : ${colors.reset} ${workflow.runName}"
        log.info "${colors.dim}  containerEngine  : ${colors.reset} ${workflow.containerEngine}"
        if (workflow.containerEngine) {
            log.info "${colors.dim}  container        : ${colors.reset} ${workflow.container}"
        }
        log.info "${colors.dim}  launchDir        : ${colors.reset} ${workflow.launchDir}"
        log.info "${colors.dim}  workDir          : ${colors.reset} ${workflow.workDir}"
        log.info "${colors.dim}  projectDir       : ${colors.reset} ${workflow.projectDir}"
        log.info "${colors.dim}  userName         : ${colors.reset} ${workflow.userName}"
        log.info "${colors.dim}  profile          : ${colors.reset} ${workflow.profile}"
        log.info "${colors.dim}  configFiles      : ${colors.reset} ${workflow.configFiles}"
        log.info NfcoreTemplate.dashedLine(params.monochrome_logs)
        log.info "${colors.green}Input/output options${colors.reset}"
        log.info "${colors.dim}  input            : ${colors.reset} ${params.input}"
        log.info "${colors.dim}  outdir           : ${colors.reset} ${params.outdir}"
        if (params.email) {
            log.info "${colors.dim}  email            : ${colors.reset} ${params.email}"
        }
        log.info NfcoreTemplate.dashedLine(params.monochrome_logs)
    }

    //
    // Dashedline
    //
    public static String dashedLine(monochrome_logs) {
        return NfcoreTemplate.dashedLine(monochrome_logs)
    }
}

// END OF SCRIPT
