process BEDTOOLS_INTERSECT_CRM {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(peaks)
    path crm_bed
    val overlap_fraction

    output:
    tuple val(meta), path("*.crm_intersected.bed"), emit: intersected
    tuple val(meta), path("*.crm_non_intersected.bed"), emit: non_intersected
    tuple val(meta), path("*.crm_annotation.bed"), emit: annotation
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def overlap_arg = overlap_fraction > 0 ? "-f ${overlap_fraction}" : ""
    """
    # Intersect peaks with CRM regions
    bedtools intersect \\
        -a $peaks \\
        -b $crm_bed \\
        -wa -wb \\
        $overlap_arg \\
        > ${prefix}.crm_intersected_raw.bed

    # Process intersected peaks and add CRM annotation
    intersect_crm_annotate.py \\
        --intersected ${prefix}.crm_intersected_raw.bed \\
        --peaks $peaks \\
        --output_intersected ${prefix}.crm_intersected.bed \\
        --output_non_intersected ${prefix}.crm_non_intersected.bed \\
        --output_annotation ${prefix}.crm_annotation.bed \\
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
    touch ${prefix}.crm_intersected.bed
    touch ${prefix}.crm_non_intersected.bed
    touch ${prefix}.crm_annotation.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

# END OF MODULE
