import nextflow.Nextflow

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
            log.info "${workflow.manifest.name} ${workflow.manifest.version}"
            System.exit(0)
        }

        // Print citation and exit on --help
        if (params.help) {
            def command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --crm_bed crm.bed --intron_bed introns.bed --gtf annotation.gtf -profile docker"
            log.info citation(workflow)
            log.info "----------------------------------------------------"
            log.info "USAGE:"
            log.info "  The typical command for running the pipeline is as follows:"
            log.info "  ${command}"
            log.info "----------------------------------------------------"
            System.exit(0)
        }

        // Check that conda channels are set-up correctly
        if (workflow.profile.contains('conda')) {
            try {
                Utils.checkCondaChannels(log)
            } catch (Exception e) {
                log.debug "Could not check conda channels: ${e.message}"
            }
        }

        // Check input has been provided
        if (!params.input) {
            Nextflow.error("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
        }
    }

    //
    // Generate methods description for MultiQC
    //
    public static String toolCitationText(params) {
        def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016),",
            "HOMER (Heinz et al. 2010),",
            "BEDTools (Quinlan and Hall 2010),",
            "MACS2 (Zhang et al. 2008)."
        ].join(' ').trim()

        return citation_text
    }

    //
    // Print parameter summary log to screen
    //
    public static void paramsSummaryLog(workflow, params, log) {
        log.info "----------------------------------------------------"
        log.info "Core Nextflow options"
        if (workflow.revision) {
            log.info "  revision         : ${workflow.revision}"
        }
        log.info "  runName          : ${workflow.runName}"
        log.info "  containerEngine  : ${workflow.containerEngine}"
        log.info "  launchDir        : ${workflow.launchDir}"
        log.info "  workDir          : ${workflow.workDir}"
        log.info "  projectDir       : ${workflow.projectDir}"
        log.info "  userName         : ${workflow.userName}"
        log.info "  profile          : ${workflow.profile}"
        log.info "  configFiles      : ${workflow.configFiles}"
        log.info "----------------------------------------------------"
        log.info "Input/output options"
        log.info "  input            : ${params.input}"
        log.info "  outdir           : ${params.outdir}"
        if (params.email) {
            log.info "  email            : ${params.email}"
        }
        log.info "----------------------------------------------------"
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</span>'}</samp></dd>\n"
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
    // Simple email function
    //
    public static void email(workflow, params, summary_params, projectDir, log, multiqc_report=[]) {
        def subject = "${workflow.manifest.name} ${workflow.success ? 'Successful' : 'FAILED'}: $workflow.runName"

        log.info "Pipeline completed: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
        log.info "Results are available in: ${params.outdir}"

        // Write simple summary
        def output_d = new File("${params.outdir}/pipeline_info/")
        if (!output_d.exists()) {
            output_d.mkdirs()
        }
        def output_hf = new File(output_d, "pipeline_report.html")
        output_hf.withWriter { w ->
            w << "<html><head><title>Pipeline Report</title></head><body>"
            w << "<h1>Pipeline Report</h1>"
            w << "<p><strong>Status:</strong> ${workflow.success ? 'SUCCESS' : 'FAILED'}</p>"
            w << "<p><strong>Run Name:</strong> ${workflow.runName}</p>"
            w << "<p><strong>Started:</strong> ${workflow.start}</p>"
            w << "<p><strong>Completed:</strong> ${workflow.complete}</p>"
            w << "<p><strong>Duration:</strong> ${workflow.duration}</p>"
            w << "</body></html>"
        }
    }
}
