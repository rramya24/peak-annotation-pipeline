process PREPARE_ANNOTATION_OUTPUT {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::python=3.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(gene_symbol_files)
    tuple val(meta2), path(consensus_peaks)

    output:
    tuple val(meta), path("*.annotation_report.html"), emit: report
    tuple val(meta), path("*.all_target_genes.txt"), emit: all_genes
    tuple val(meta), path("*.peak_annotation_summary.txt"), emit: summary
    tuple val(meta), path("*.mqc.tsv"), emit: multiqc_files
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Prepare final annotation output
    prepare_annotation_output.py \\
        --gene_symbol_files ${gene_symbol_files.join(' ')} \\
        --consensus_peaks $consensus_peaks \\
        --output_prefix $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annotation_report.html
    touch ${prefix}.all_target_genes.txt
    touch ${prefix}.peak_annotation_summary.txt
    touch ${prefix}.mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


