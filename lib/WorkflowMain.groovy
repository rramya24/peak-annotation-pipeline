class WorkflowMain {
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name ?: 'this pipeline'} for your analysis please cite:\n\n* The pipeline\n  https://github.com/rramya24/peak-annotation-pipeline"
    }
    public static void initialise(workflow, params, log) {
        if (params.version) {
            log.info "${workflow.manifest.name ?: 'peak-annotation-pipeline'} ${workflow.manifest.version ?: '1.0.0'}"
            System.exit(0)
        }
        if (params.help) {
            log.info citation(workflow)
            System.exit(0)
        }
        if (workflow.profile.contains('conda')) {
            try {
                Utils.checkCondaChannels(log)
            } catch (Exception e) {
                log.debug "Could not check conda channels: ${e.message}"
            }
        }
        if (!params.input) {
            throw new Exception("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
        }
    }
    public static String workflowSummary(workflow, params) {
        return "Pipeline: ${workflow.manifest.name ?: 'peak-annotation-pipeline'}\nInput: ${params.input ?: 'Not specified'}"
    }
    public static String toolCitationText(params) {
        return "Tools used: MultiQC, HOMER, BEDTools, MACS2."
    }
    public static void paramsSummaryLog(workflow, params, log) {
        log.info "Pipeline: ${workflow.manifest.name ?: 'peak-annotation-pipeline'}"
        log.info "Input: ${params.input ?: 'Not specified'}"
    }
}
