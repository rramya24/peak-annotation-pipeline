process BEDTOOLS_INTERSECT_INTRON {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/mulled-v2-2f48cc59b03027e31ead6d383fe1b8057785dd24:5d182f583f4696f4c4d9f3be93052811b383341f-0' :
    'biocontainers/mulled-v2-2f48cc59b03027e31ead6d383fe1b8057785dd24:5d182f583f4696f4c4d9f3be93052811b383341f-0' }"

    input:
    tuple val(meta), path(peaks)
    path intron_bed
    path gtf
    val overlap_fraction

    output:
    tuple val(meta), path("*.intron_intersected.bed"), emit: intersected
    tuple val(meta), path("*.intron_non_intersected.bed"), emit: non_intersected
    tuple val(meta), path("*.intron_annotation.tsv"), emit: annotation
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gtf_param = gtf ? "--gtf $gtf" : ""
    """
    # Single Python script does everything
    $projectDir/bin/intersect_intron_annotate.py \\
        --peaks $peaks \\
        --introns $intron_bed \\
        $gtf_param \\
        --intersected ${prefix}.intron_intersected.bed \\
        --annotated ${prefix}.intron_annotation.tsv \\
        --non-intersected ${prefix}.intron_non_intersected.bed \\
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
    touch ${prefix}.intron_intersected.bed
    touch ${prefix}.intron_non_intersected.bed
    touch ${prefix}.intron_annotation.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
