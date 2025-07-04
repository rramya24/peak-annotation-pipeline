//
// Consensus peak calling subworkflow
// Groups replicates by sample name and processes each sample separately
//

include { CONSENSUS_PEAKS } from '../../modules/local/consensus_peaks/main'

workflow MACS2_CONSENSUS {
    take:
    peaks           // channel: [ val(meta), path(peaks) ]
    min_reps        // integer: minimum number of replicates

    main:
    ch_versions = Channel.empty()

    //
    // Group peaks by sample name and collect all replicates for each sample
    //
    peaks
        .map { meta, peak_file ->
            // Create a tuple with sample name and peak file info
            def sample_name = meta.id
            def replicate = meta.replicate
            return [ sample_name, [ meta, peak_file ] ]
        }
        .groupTuple()
        .map { sample_name, meta_peak_list ->
            // Create new meta for the sample (without replicate info)
            def sample_meta = [:]
            sample_meta.id = sample_name
            sample_meta.single_end = true

            // Extract peak files from the grouped data
            def peak_files = meta_peak_list.collect { it[1] }

            return [ sample_meta, peak_files ]
        }
        .set { ch_grouped_peaks }

    //
    // Call consensus peaks for each sample
    //
    CONSENSUS_PEAKS (
        ch_grouped_peaks,
        min_reps
    )
    ch_versions = ch_versions.mix(CONSENSUS_PEAKS.out.versions)

    emit:
    consensus_peaks = CONSENSUS_PEAKS.out.consensus_peaks  // channel: [ val(meta), path(consensus_peaks) ]
    log            = CONSENSUS_PEAKS.out.log              // channel: [ val(meta), path(log) ]
    versions       = ch_versions                          // channel: [ versions.yml ]
}

# END OF SUBWORKFLOW
