process CONVERT_GENEID_TO_SYMBOL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--h7d875b9_6' :
        'biocontainers/bioawk:1.0--h7d875b9_6' }"

    input:
    tuple val(meta), path(annotation_file)
    path gtf

    output:
    tuple val(meta), path("*.converted.txt"), emit: converted_annotation
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create a mapping from gene_id to gene_name/symbol from GTF
    bioawk -c gff '
        \$feature=="gene" && \$attribute ~ /gene_id/ && \$attribute ~ /gene_name/ {
            match(\$attribute, /gene_id "([^"]+)"/, gene_id)
            match(\$attribute, /gene_name "([^"]+)"/, gene_name)
            if (gene_id[1] && gene_name[1]) {
                print gene_id[1] "\\t" gene_name[1]
            }
        }
    ' $gtf > gene_id_to_symbol.map

    # Convert gene IDs to symbols in annotation file
    # Assuming the annotation file has gene IDs in a specific column
    # This will need to be adjusted based on your actual file format
    awk 'BEGIN{OFS="\\t"}
        NR==FNR {map[\$1]=\$2; next}
        {
            # Replace gene IDs with symbols where found
            for (i=1; i<=NF; i++) {
                if (\$i in map) {
                    \$i = map[\$i]
                }
            }
            print
        }' gene_id_to_symbol.map $annotation_file > ${prefix}.converted.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioawk: \$(bioawk --version 2>&1 | head -n1 | sed 's/^.*bioawk //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.converted.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioawk: \$(bioawk --version 2>&1 | head -n1 | sed 's/^.*bioawk //; s/ .*\$//')
    END_VERSIONS
    """
}
