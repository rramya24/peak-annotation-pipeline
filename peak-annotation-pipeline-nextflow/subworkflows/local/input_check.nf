//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_peaks_channel(it) }
        .set { peaks }

    emit:
    peaks                                     // channel: [ val(meta), [ peaks ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ peaks ] ]
def create_peaks_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id                = row.sample
    meta.single_end        = true
    meta.replicate         = row.replicate.toInteger()

    // add path(s) of the peak file(s) to the meta map
    def peaks_meta = []
    if (!file(row.peaks).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Peak file does not exist!\n${row.peaks}"
    }

    // Create a processed peaks file with proper identifiers
    def processed_peaks = add_peak_identifiers(file(row.peaks), row.sample, row.replicate)
    peaks_meta = [ meta, processed_//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_peaks_channel(it) }
        .set { peaks }

    emit:
    peaks                                     // channel: [ val(meta), [ peaks ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ peaks ] ]
def create_peaks_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id                = row.sample
    meta.single_end        = true
    meta.replicate         = row.replicate.toInteger()

    // add path(s) of the peak file(s) to the meta map
    def peaks_meta = []
    if (!file(row.peaks).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Peak file does not exist!\n${row.peaks}"
    }

    peaks_meta = [ meta, file(row.peaks) ]
    return peaks_meta
}

# END OF SUBWORKFLOW
