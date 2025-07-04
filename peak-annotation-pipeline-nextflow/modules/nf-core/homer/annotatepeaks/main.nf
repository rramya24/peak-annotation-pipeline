process HOMER_ANNOTATEPEAKS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::homer=4.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl5321h9a82719_6' :
        'quay.io/biocontainers/homer:4.11--pl5321h9a82719_6' }"

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
    """
    # This matches the exact HOMER module from nf-core/chipseq
    annotatePeaks.pl \\
        $peaks \\
        $genome \\
        -gtf $gtf \\
        $distance_arg \\
        $args \\
        > ${prefix}.annotatePeaks.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: \$(echo \$(homer2 --version 2>&1) | sed 's/^.*homer2 //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annotatePeaks.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: \$(echo \$(homer2 --version 2>&1) | sed 's/^.*homer2 //; s/Using.*\$//')
    END_VERSIONS
    """
}

# END OF MODULE
