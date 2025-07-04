process BEDTOOLS_INTERSECT_INTRON {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(peaks)
    path intron_bed
    val overlap_fraction

    output:
    tuple val(meta), path("*.intron_intersected.bed"), emit: intersected
    tuple val(meta), path("*.intron_non_intersected.bed"), emit: non_intersected
    tuple val(meta), path("*.intron_annotation.bed"), emit: annotation
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def overlap_arg = overlap_fraction > 0 ? "-f ${overlap_fraction}" : ""
    """
    # Intersect peaks with intron regions
    bedtools intersect \\
        -a $peaks \\
        -b $intron_bed \\
        -wa -wb \\
        $overlap_arg \\
        > ${prefix}.intron_intersected_raw.bed

    # Process intersected peaks and add intron annotation
    intersect_intron_annotate.py \\
        --intersected ${prefix}.intron_intersected_raw.bed \\
        --peaks $peaks \\
        --output_intersected ${prefix}.intron_intersected.bed \\
        --output_non_intersected ${prefix}.intron_non_intersected.bed \\
        --output_annotation ${prefix}.intron_annotation.bed \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.intron_intersected.bed
    touch ${prefix}.intron_non_intersected.bed
    touch ${prefix}.intron_annotation.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

# END OF MODULE
