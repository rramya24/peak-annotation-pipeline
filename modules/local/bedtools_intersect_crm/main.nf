process BEDTOOLS_INTERSECT_CRM {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/mulled-v2-2f48cc59b03027e31ead6d383fe1b8057785dd24:5d182f583f4696f4c4d9f3be93052811b383341f-0' :
    'biocontainers/mulled-v2-2f48cc59b03027e31ead6d383fe1b8057785dd24:5d182f583f4696f4c4d9f3be93052811b383341f-0' }"

    input:
    tuple val(meta), path(peaks)
    path crm_bed
    path gtf
    val overlap_fraction

    output:
    tuple val(meta), path("*.crm_intersected.bed"), emit: intersected
    tuple val(meta), path("*.crm_non_intersected.bed"), emit: non_intersected
    tuple val(meta), path("*.crm_annotation.tsv"), emit: annotation
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gtf_param = gtf ? "--gtf $gtf" : ""
    """
    # Single Python script does everything
    $projectDir/bin/intersect_crm_annotate.py \\
        --peaks $peaks \\
        --crm $crm_bed \\
        $gtf_param \\
        --intersected ${prefix}.crm_intersected.bed \\
        --annotated ${prefix}.crm_annotation.tsv \\
        --non-intersected ${prefix}.crm_non_intersected.bed \\
        --overlap-fraction $overlap_fraction

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
    touch ${prefix}.crm_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
