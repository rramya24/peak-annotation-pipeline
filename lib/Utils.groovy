//
// This file holds several utility functions used within the nf-core pipeline template.
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class Utils {

    //
    // When running with -profile conda, warn if channels have not been set-up appropriately
    //
    public static void checkCondaChannels(log) {
        Nextflow.script('bash', ['-c', 'conda config --show channels']).out.splitText().each { line ->
            if (line.startsWith('channels:')) {
                def channels = []
                return
            }
            if (line.startsWith('  - ') && line.contains('conda-forge')) {
                channels << 'conda-forge'
            }
            if (line.startsWith('  - ') && line.contains('bioconda')) {
                channels << 'bioconda'
            }
        }
        if (!channels.contains('conda-forge')) {
            log.warn "conda-forge channel not found. Add with: conda config --add channels conda-forge"
        }
        if (!channels.contains('bioconda')) {
            log.warn "bioconda channel not found. Add with: conda config --add channels bioconda"
        }
    }

    //
    // Check if a row has the expected number of columns
    //
    public static Boolean checkNumberOfItemsInRow(row, number_of_columns, log) {
        if (row.size() != number_of_columns) {
            log.error("Invalid number of columns (minimum = ${number_of_columns}) in row: ${row}")
            return false
        }
        return true
    }

    //
    // Return file extension
    //
    public static String getFileExtension(filename) {
        def ext = ""
        if (filename.indexOf(".") != -1) {
            ext = filename.substring(filename.lastIndexOf("."))
        }
        return ext
    }

    //
    // Remove file extension from filename
    //
    public static String removeFileExtension(filename) {
        def ext = getFileExtension(filename)
        if (ext) {
            return filename.substring(0, filename.lastIndexOf(ext))
        }
        return filename
    }

    //
    // Check if file exists
    //
    public static Boolean checkFileExists(file_path, log) {
        if (!file(file_path).exists()) {
            log.error("File does not exist: ${file_path}")
            return false
        }
        return true
    }

    //
    // Check if file is readable
    //
    public static Boolean checkFileReadable(file_path, log) {
        if (!file(file_path).canRead()) {
            log.error("File is not readable: ${file_path}")
            return false
        }
        return true
    }

    //
    // Check if directory exists
    //
    public static Boolean checkDirectoryExists(dir_path, log) {
        if (!file(dir_path).exists()) {
            log.error("Directory does not exist: ${dir_path}")
            return false
        }
        if (!file(dir_path).isDirectory()) {
            log.error("Path is not a directory: ${dir_path}")
            return false
        }
        return true
    }

    //
    // Function to generate a UUID
    //
    public static String generateUUID() {
        return UUID.randomUUID().toString()
    }

    //
    // Function to get basename of file
    //
    public static String getBasename(filename) {
        def file_obj = new File(filename)
        return file_obj.getName()
    }

    //
    // Function to get parent directory of file
    //
    public static String getParentDir(filename) {
        def file_obj = new File(filename)
        return file_obj.getParent() ?: '.'
    }

    //
    // Function to create directory if it doesn't exist
    //
    public static void createDirectory(dir_path, log) {
        def dir = new File(dir_path)
        if (!dir.exists()) {
            if (dir.mkdirs()) {
                log.info("Created directory: ${dir_path}")
            } else {
                log.error("Failed to create directory: ${dir_path}")
            }
        }
    }

    //
    // Function to format duration
    //
    public static String formatDuration(duration) {
        def seconds = duration.toSeconds()
        def hours = Math.floor(seconds / 3600)
        def minutes = Math.floor((seconds % 3600) / 60)
        def secs = seconds % 60

        if (hours > 0) {
            return String.format("%02d:%02d:%02d", hours, minutes, secs)
        } else {
            return String.format("%02d:%02d", minutes, secs)
        }
    }

    //
    // Function to format memory
    //
    public static String formatMemory(memory) {
        if (memory.toString().contains('GB')) {
            return memory.toString()
        } else if (memory.toString().contains('MB')) {
            return memory.toString()
        } else if (memory.toString().contains('KB')) {
            return memory.toString()
        } else {
            // Assume bytes
            def mb = memory / (1024 * 1024)
            if (mb >= 1024) {
                def gb = mb / 1024
                return String.format("%.1f GB", gb)
            } else {
                return String.format("%.1f MB", mb)
            }
        }
    }

    //
    // Function to get current timestamp
    //
    public static String getCurrentTimestamp() {
        return new Date().format("yyyy-MM-dd HH:mm:ss")
    }

    //
    // Function to safely get parameter value
    //
    public static Object getParamValue(params, param_name, default_value = null) {
        if (params.containsKey(param_name)) {
            return params[param_name]
        }
        return default_value
    }

    //
    // Function to validate email address
    //
    public static Boolean isValidEmail(email) {
        def pattern = /^[A-Za-z0-9+_.-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}$/
        return email.matches(pattern)
    }

    //
    // Function to clean filename for safe usage
    //
    public static String cleanFilename(filename) {
        return filename.replaceAll(/[^a-zA-Z0-9._-]/, '_')
    }

    //
    // Function to join paths safely
    //
    public static String joinPaths(String... paths) {
        return paths.join(File.separator)
    }

    //
    // Function to check if string is numeric
    //
    public static Boolean isNumeric(str) {
        try {
            Double.parseDouble(str)
            return true
        } catch (NumberFormatException e) {
            return false
        }
    }

    //
    // Function to truncate string to specified length
    //
    public static String truncateString(str, maxLength) {
        if (str.length() <= maxLength) {
            return str
        }
        return str.substring(0, maxLength - 3) + "..."
    }

    //
    // Function to capitalize first letter
    //
    public static String capitalize(str) {
        if (!str || str.length() == 0) {
            return str
        }
        return str[0].toUpperCase() + str.substring(1)
    }

    //
    // Function to convert bytes to human readable format
    //
    public static String bytesToHuman(bytes) {
        def units = ['B', 'KB', 'MB', 'GB', 'TB']
        def unitIndex = 0
        def size = bytes.toDouble()

        while (size >= 1024 && unitIndex < units.size() - 1) {
            size /= 1024
            unitIndex++
        }

        return String.format("%.2f %s", size, units[unitIndex])
    }

    //
    // Function to merge maps
    //
    public static Map mergeMaps(Map... maps) {
        def result = [:]
        for (map in maps) {
            if (map) {
                result.putAll(map)
            }
        }
        return result
    }

    //
    // Function to check if process should be skipped
    //
    public static Boolean skipProcess(params, process_name) {
        def skip_param = "skip_${process_name}"
        return params.containsKey(skip_param) && params[skip_param]
    }

    //
    // Function to get process resources
    //
    public static Map getProcessResources(params, process_name) {
        def resources = [:]

        // CPU
        def cpu_param = "${process_name}_cpu"
        if (params.containsKey(cpu_param)) {
            resources.cpus = params[cpu_param]
        }

        // Memory
        def memory_param = "${process_name}_memory"
        if (params.containsKey(memory_param)) {
            resources.memory = params[memory_param]
        }

        // Time
        def time_param = "${process_name}_time"
        if (params.containsKey(time_param)) {
            resources.time = params[time_param]
        }

        return resources
    }

    //
    // Function to validate file extension
    //
    public static Boolean hasValidExtension(filename, validExtensions) {
        def ext = getFileExtension(filename).toLowerCase()
        return validExtensions.any { it.toLowerCase() == ext }
    }
}

// END OF SCRIPT
