process CONSENSUS_PEAKS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(peaks)
    val min_reps

    output:
    tuple val(meta), path("*.consensus_peaks.bed"), emit: consensus_peaks
    tuple val(meta), path("*.consensus_peaks.boolean.intersect.plot.pdf"), emit: pdf, optional: true
    tuple val(meta), path("*.consensus_peaks.boolean.intersect.txt"), emit: txt, optional: true
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def peak_files = peaks.collect{it.toString()}.join(' ')
    """
    # This is the exact same logic from nf-core/chipseq v2.0
    sort -k1,1V -k2,2n -k3,3n ${peak_files} | mergeBed -i stdin -d 150 -c 4 -o count_distinct > ${prefix}.consensus_peaks.txt

    consensus_peaks.py \\
        --peak_files ${peak_files} \\
        --min_reps $min_reps \\
        --prefix $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.consensus_peaks.bed
    touch ${prefix}.consensus_peaks.boolean.intersect.plot.pdf
    touch ${prefix}.consensus_peaks.boolean.intersect.txt
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}


