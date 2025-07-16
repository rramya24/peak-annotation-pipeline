process EXPAND_TARGETS_LNCRNA {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::python=3.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(target_genes)
    path lncrna_mirna_mapping

    output:
    tuple val(meta), path("*.expanded_targets.txt"), emit: expanded_targets
    tuple val(meta), path("*.expansion_log.txt"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Expand target genes with lncRNA-miRNA relationships
    expand_targets_with_lncrna_mirna.py \\
        --target_genes $target_genes \\
        --lncrna_mirna_mapping $lncrna_mirna_mapping \\
        --output ${prefix}.expanded_targets.txt \\
        --log ${prefix}.expansion_log.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.expanded_targets.txt
    touch ${prefix}.expansion_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


