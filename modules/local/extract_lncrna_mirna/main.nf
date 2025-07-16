process EXTRACT_LNCRNA_MIRNA {
    tag "extract_lncrna_mirna"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::python=3.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    path gtf

    output:
    path "lncrna_mirna_mapping.txt", emit: mapping
    path "lncrna_mirna_extraction.log", emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Extract lncRNA-miRNA relationships from GTF file
    extract_lncrna_mirna_relationships.py \\
        --gtf $gtf \\
        --output lncrna_mirna_mapping.txt \\
        --log lncrna_mirna_extraction.log \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch lncrna_mirna_mapping.txt
    touch lncrna_mirna_extraction.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


