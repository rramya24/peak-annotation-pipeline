//
// This file holds several functions used for validation of parameters by nf-core pipelines
//

import nextflow.Nextflow
import org.yaml.snakeyaml.Yaml
import groovy.json.JsonSlurper

class NfcoreSchema {

    //
    // Resolve Schema path relative to main workflow directory
    //
    public static String getSchemaPath(workflow, schema_filename='nextflow_schema.json') {
        return "${workflow.projectDir}/${schema_filename}"
    }

    //
    // Function to loop over all parameters defined in schema and check
    // whether the given parameters adhere to the specifications
    //
    public static void validateParameters(workflow, params, log, schema_filename='nextflow_schema.json') {
        def schema_path = getSchemaPath(workflow, schema_filename)
        def schema_file = new File(schema_path)
        if (!schema_file.exists()) {
            log.warn "Schema file not found: ${schema_path}"
            return
        }

        def slurper = new JsonSlurper()
        def schema = slurper.parse(schema_file)

        // Validate top-level parameters
        def schema_params = schema.get('definitions', [:])
        for (group in schema_params.keySet()) {
            def group_params = schema_params[group].get('properties', [:])
            for (param in group_params.keySet()) {
                def param_spec = group_params[param]
                def param_value = params[param]

                // Skip if parameter not provided
                if (param_value == null) continue

                // Check parameter type
                if (param_spec.containsKey('type')) {
                    validateParameterType(param, param_value, param_spec.type, log)
                }

                // Check parameter format
                if (param_spec.containsKey('format')) {
                    validateParameterFormat(param, param_value, param_spec.format, log)
                }

                // Check parameter pattern
                if (param_spec.containsKey('pattern')) {
                    validateParameterPattern(param, param_value, param_spec.pattern, log)
                }

                // Check parameter minimum/maximum
                if (param_spec.containsKey('minimum') && param_value instanceof Number) {
                    if (param_value < param_spec.minimum) {
                        log.error "Parameter '${param}' value ${param_value} is below minimum ${param_spec.minimum}"
                    }
                }
                if (param_spec.containsKey('maximum') && param_value instanceof Number) {
                    if (param_value > param_spec.maximum) {
                        log.error "Parameter '${param}' value ${param_value} is above maximum ${param_spec.maximum}"
                    }
                }

                // Check parameter enum
                if (param_spec.containsKey('enum')) {
                    if (!param_spec.enum.contains(param_value)) {
                        log.error "Parameter '${param}' value '${param_value}' is not in allowed values: ${param_spec.enum}"
                    }
                }
            }
        }
    }

    //
    // Validate parameter type
    //
    public static void validateParameterType(param, value, expected_type, log) {
        def actual_type = value.getClass().getSimpleName().toLowerCase()

        // Map Java types to JSON schema types
        def type_map = [
            'string': 'string',
            'integer': 'integer',
            'biginteger': 'integer',
            'double': 'number',
            'bigdecimal': 'number',
            'float': 'number',
            'boolean': 'boolean',
            'arraylist': 'array',
            'linkedhashmap': 'object'
        ]

        def mapped_type = type_map[actual_type] ?: actual_type

        if (mapped_type != expected_type) {
            log.error "Parameter '${param}' expected type '${expected_type}' but got '${mapped_type}'"
        }
    }

    //
    // Validate parameter format
    //
    public static void validateParameterFormat(param, value, format, log) {
        switch (format) {
            case 'file-path':
                if (!value.toString().startsWith('/') && !value.toString().startsWith('./') && !value.toString().startsWith('../')) {
                    log.warn "Parameter '${param}' should be a file path but got: ${value}"
                }
                break
            case 'directory-path':
                if (!value.toString().startsWith('/') && !value.toString().startsWith('./') && !value.toString().startsWith('../')) {
                    log.warn "Parameter '${param}' should be a directory path but got: ${value}"
                }
                break
            case 'path':
                if (!value.toString().startsWith('/') && !value.toString().startsWith('./') && !value.toString().startsWith('../')) {
                    log.warn "Parameter '${param}' should be a path but got: ${value}"
                }
                break
            case 'email':
                if (!value.toString().contains('@')) {
                    log.error "Parameter '${param}' should be an email address but got: ${value}"
                }
                break
            default:
                // No specific validation for other formats
                break
        }
    }

    //
    // Validate parameter pattern
    //
    public static void validateParameterPattern(param, value, pattern, log) {
        try {
            if (!value.toString().matches(pattern)) {
                log.error "Parameter '${param}' value '${value}' does not match pattern '${pattern}'"
            }
        } catch (Exception e) {
            log.warn "Could not validate pattern for parameter '${param}': ${e.message}"
        }
    }

    //
    // Function to retrieve the nextflow schema in JSON format
    //
    public static Map getSchemaMap(workflow, schema_filename='nextflow_schema.json') {
        def schema_path = getSchemaPath(workflow, schema_filename)
        def schema_file = new File(schema_path)
        if (!schema_file.exists()) {
            return [:]
        }

        def slurper = new JsonSlurper()
        return slurper.parse(schema_file)
    }

    //
    // Function to collect the pipeline parameters with their defaults
    //
    public static LinkedHashMap paramsSummaryMap(workflow, params, schema_filename='nextflow_schema.json') {
        // Get the schema
        def schema = getSchemaMap(workflow, schema_filename)

        // Collect parameters grouped by schema sections
        def summary_params = [:]

        // Core Nextflow options
        summary_params['Core Nextflow options'] = [:]
        summary_params['Core Nextflow options']['revision'] = workflow.revision ?: 'None'
        summary_params['Core Nextflow options']['runName'] = workflow.runName
        summary_params['Core Nextflow options']['containerEngine'] = workflow.containerEngine ?: 'None'
        summary_params['Core Nextflow options']['container'] = workflow.container ?: 'None'
        summary_params['Core Nextflow options']['launchDir'] = workflow.launchDir
        summary_params['Core Nextflow options']['workDir'] = workflow.workDir
        summary_params['Core Nextflow options']['projectDir'] = workflow.projectDir
        summary_params['Core Nextflow options']['userName'] = workflow.userName
        summary_params['Core Nextflow options']['profile'] = workflow.profile
        summary_params['Core Nextflow options']['configFiles'] = workflow.configFiles.join(', ')

        // Get parameters from schema
        def schema_params = schema.get('definitions', [:])
        for (group in schema_params.keySet()) {
            def group_name = schema_params[group].get('title', group)
            def group_description = schema_params[group].get('description', '')
            def group_params = schema_params[group].get('properties', [:])

            if (group_params.size() > 0) {
                summary_params[group_name] = [:]
                for (param in group_params.keySet()) {
                    def param_value = params[param]
                    def param_default = group_params[param].get('default', null)

                    // Use provided value or default
                    if (param_value != null) {
                        summary_params[group_name][param] = param_value
                    } else if (param_default != null) {
                        summary_params[group_name][param] = param_default
                    } else {
                        summary_params[group_name][param] = null
                    }
                }
            }
        }

        return summary_params
    }

    //
    // Function to print a parameter summary map
    //
    public static String paramsSummaryLog(workflow, params, log, schema_filename='nextflow_schema.json') {
        Map colors = NfcoreTemplate.logColours(params.monochrome_logs)
        String output = ''

        def summary_params = paramsSummaryMap(workflow, params, schema_filename)

        for (group in summary_params.keySet()) {
            def group_params = summary_params[group]
            if (group_params.size() > 0) {
                output += "${colors.green}${group}${colors.reset}\n"
                for (param in group_params.keySet()) {
                    def param_value = group_params[param]
                    def param_display = param_value ?: "${colors.dim}None${colors.reset}"
                    if (param_value instanceof Boolean) {
                        param_display = param_value ? "${colors.green}true${colors.reset}" : "${colors.red}false${colors.reset}"
                    }
                    output += "${colors.dim}  ${param.padRight(25)}: ${colors.reset}${param_display}\n"
                }
                output += "\n"
            }
        }

        return output
    }

    //
    // Function to validate a sample sheet
    //
    public static void validateSampleSheet(samplesheet_file, log) {
        if (!samplesheet_file.exists()) {
            log.error "Sample sheet file not found: ${samplesheet_file}"
            return
        }

        // Basic CSV validation
        def lines = samplesheet_file.readLines()
        if (lines.size() < 2) {
            log.error "Sample sheet must contain at least a header and one data row"
            return
        }

        // Check header
        def header = lines[0].split(',').collect { it.trim() }
        def required_columns = ['sample', 'replicate', 'peak_file']

        for (column in required_columns) {
            if (!header.contains(column)) {
                log.error "Sample sheet missing required column: ${column}"
            }
        }

        // Validate data rows
        for (i in 1..<lines.size()) {
            def row = lines[i].split(',').collect { it.trim() }
            if (row.size() != header.size()) {
                log.error "Sample sheet row ${i + 1} has ${row.size()} columns but header has ${header.size()}"
            }
        }

        log.info "Sample sheet validation completed"
    }

    //
    // Function to get default parameters from schema
    //
    public static Map getDefaultParams(workflow, schema_filename='nextflow_schema.json') {
        def schema = getSchemaMap(workflow, schema_filename)
        def default_params = [:]

        def schema_params = schema.get('definitions', [:])
        for (group in schema_params.keySet()) {
            def group_params = schema_params[group].get('properties', [:])
            for (param in group_params.keySet()) {
                def param_spec = group_params[param]
                if (param_spec.containsKey('default')) {
                    default_params[param] = param_spec.default
                }
            }
        }

        return default_params
    }

    //
    // Function to create a parameters file
    //
    public static void createParamsFile(workflow, params, outfile, schema_filename='nextflow_schema.json') {
        def schema = getSchemaMap(workflow, schema_filename)
        def yaml = new Yaml()

        // Create YAML content
        def yaml_content = [:]

        def schema_params = schema.get('definitions', [:])
        for (group in schema_params.keySet()) {
            def group_params = schema_params[group].get('properties', [:])
            for (param in group_params.keySet()) {
                def param_value = params[param]
                if (param_value != null) {
                    yaml_content[param] = param_value
                }
            }
        }

        // Write to file
        def output_file = new File(outfile)
        output_file.withWriter { writer ->
            yaml.dump(yaml_content, writer)
        }
    }
}

// END OF SCRIPT
