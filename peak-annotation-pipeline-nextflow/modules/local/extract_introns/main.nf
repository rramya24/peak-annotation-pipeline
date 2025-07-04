process EXTRACT_INTRONS {
    tag "extract_introns"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::python=3.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    path gtf

    output:
    path "first_introns.bed", emit: introns
    path "intron_extraction.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Extract first introns from GTF file
    extract_first_introns.py \\
        --gtf $gtf \\
        --output first_introns.bed \\
        --log intron_extraction.log \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch first_introns.bed
    touch intron_extraction.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

# END OF MODULE
