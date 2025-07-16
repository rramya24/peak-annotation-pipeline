//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    // Check samplesheet format and validate contents
    SAMPLESHEET_CHECK ( samplesheet )

    // Extract valid sample information and create channels
    ch_samples = SAMPLESHEET_CHECK.out.csv
        .splitCsv(header:true, sep:',')
        .map { row -> create_peaks_channel(row) }

    // Group by sample ID for consensus calling (multi-replicate support)
    ch_samples_grouped = ch_samples
        .groupTuple(by: 0)
        .map { meta, peaks_files ->
            // Create grouped meta with replicate information
            def grouped_meta = meta.clone()
            grouped_meta.replicates = peaks_files.size()
            grouped_meta.replicate_ids = peaks_files.collect { it.name }

            return [grouped_meta, peaks_files]
        }

    // Individual samples for processes that need single files
    ch_individual_samples = ch_samples
        .map { meta, peaks_file ->
            // Add replicate-specific ID
            def rep_meta = meta.clone()
            rep_meta.replicate_id = "${meta.id}_rep${meta.replicate}"

            return [rep_meta, peaks_file]
        }

    emit:
    samples = ch_samples_grouped        // channel: [ meta, [peaks_files] ]
    individual = ch_individual_samples  // channel: [ meta, peaks_file ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to create peaks channel with validation
def create_peaks_channel(LinkedHashMap row) {
    // Validate required columns
    if (!row.sample) {
        log.error "Sample name is required in samplesheet"
        System.exit(1)
    }
    if (!row.replicate) {
        log.error "Replicate number is required in samplesheet"
        System.exit(1)
    }
    if (!row.peaks) {
        log.error "Peaks file path is required in samplesheet"
        System.exit(1)
    }

    // Validate file existence
    def peaks_file = file(row.peaks)
    if (!peaks_file.exists()) {
        log.error "Peaks file does not exist: ${row.peaks}"
        System.exit(1)
    }

    // Create meta map
    def meta = [:]
    meta.id = row.sample
    meta.replicate = row.replicate.toInteger()
    meta.single_end = true // Peak files are single-end by nature

    // Return tuple with meta and peaks file
    return [meta, peaks_file]
}
