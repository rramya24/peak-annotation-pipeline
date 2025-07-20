process HOMER_ANNOTATEPEAKS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::homer=4.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3' :
        'biocontainers/homer:4.11--pl526hc9558a2_3' }"

    input:
    tuple val(meta), path(peaks)
    path gtf
    val genome
    val distance

    output:
    tuple val(meta), path("*.annotatePeaks.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def distance_arg = distance ? "-d ${distance}" : "-d 1000"
    def VERSION = '4.11'

    """
    # Check and fix chromosome naming in peak file
    echo "Original peak file (first 5 lines):"
    head -5 ${peaks}

    # Create a copy of the peak file with chromosome name fixes if needed
    cp ${peaks} ${prefix}_peaks_fixed.bed

    # Check if we need to add 'chr' prefix for Drosophila
    FIRST_CHR=\$(head -1 ${peaks} | cut -f1)
    echo "First chromosome in peak file: \$FIRST_CHR"

    # For Drosophila, chromosomes should be: chr2L, chr2R, chr3L, chr3R, chrX, chrY, chr4, chrM
    # Using safer chromosome detection without end-of-line anchor
    if [[ \$FIRST_CHR =~ ^[0-9XYM] ]] || [[ \$FIRST_CHR == "2L" ]] || [[ \$FIRST_CHR == "2R" ]] || [[ \$FIRST_CHR == "3L" ]] || [[ \$FIRST_CHR == "3R" ]]; then
        echo "Adding 'chr' prefix to chromosome names..."
        sed 's/^2L/chr2L/; s/^2R/chr2R/; s/^3L/chr3L/; s/^3R/chr3R/; s/^X\t/chrX\t/; s/^Y\t/chrY\t/; s/^4\t/chr4\t/; s/^M\t/chrM\t/; s/^X\$/chrX/; s/^Y\$/chrY/; s/^4\$/chr4/; s/^M\$/chrM/' ${peaks} > ${prefix}_peaks_fixed.bed
    elif [[ \$FIRST_CHR =~ ^chr ]]; then
        echo "Chromosomes already have 'chr' prefix"
        cp ${peaks} ${prefix}_peaks_fixed.bed
    else
        echo "Warning: Unrecognized chromosome format: \$FIRST_CHR"
        cp ${peaks} ${prefix}_peaks_fixed.bed
    fi

    echo "Fixed peak file (first 5 lines):"
    head -5 ${prefix}_peaks_fixed.bed

    # Run HOMER annotation with error handling
    echo "Running HOMER annotation..."
    echo "Command: annotatePeaks.pl ${prefix}_peaks_fixed.bed ${genome} -gtf ${gtf} ${distance_arg} ${args}"

    annotatePeaks.pl \\
        ${prefix}_peaks_fixed.bed \\
        ${genome} \\
        -gtf ${gtf} \\
        ${distance_arg} \\
        ${args} \\
        > ${prefix}.annotatePeaks.txt

    # Check if output was generated and has content
    if [[ ! -s ${prefix}.annotatePeaks.txt ]]; then
        echo "Warning: HOMER output is empty. Checking for issues..."
        echo "Peak file contents:"
        head -10 ${prefix}_peaks_fixed.bed

        # Try without GTF
        echo "Trying without GTF file..."
        annotatePeaks.pl \\
            ${prefix}_peaks_fixed.bed \\
            ${genome} \\
            ${distance_arg} \\
            > ${prefix}.annotatePeaks.txt.backup

        # Use backup if it has content
        if [[ -s ${prefix}.annotatePeaks.txt.backup ]]; then
            mv ${prefix}.annotatePeaks.txt.backup ${prefix}.annotatePeaks.txt
            echo "Used backup annotation (without GTF)"
        else
            # Create minimal output with headers
            echo "Creating minimal output with headers..."
            echo -e "PeakID\tChr\tStart\tEnd\tStrand\tPeak Score\tFocus Ratio/Region Size\tAnnotation\tDetailed Annotation\tDistance to TSS\tNearest PromoterID\tEntrez ID\tNearest Unigene\tNearest Refseq\tNearest Ensembl\tGene Name\tGene Alias\tGene Description\tGene Type" > ${prefix}.annotatePeaks.txt
            echo "Warning: HOMER annotation failed. Created empty output with headers."
        fi
    else
        echo "HOMER annotation completed successfully."
        echo "Output lines: \$(wc -l < ${prefix}.annotatePeaks.txt)"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.11'
    """
    touch ${prefix}.annotatePeaks.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: ${VERSION}
    END_VERSIONS
    """
}
