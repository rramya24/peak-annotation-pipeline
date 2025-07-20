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
    """
    # Process peaks with CRM regions using custom script
    # Script filters out 'Unspecified' annotations internally
    $projectDir/bin/intersect_crm_annotate.py \\
        --peaks $peaks \\
        --crm $crm_bed \\
        --intersected ${prefix}.crm_intersected.bed \\
        --annotated ${prefix}.crm_annotation.bed \\
        --overlap-fraction $overlap_fraction \\

    # Create non-intersected file: peaks that intersected ONLY unspecified CRMs
    # OR peaks that didn't intersect any CRMs at all
    # Script handles this logic internally by reading the intersected file

    # Get peaks that intersected with specified CRMs
    if [ -s "${prefix}.crm_intersected.bed" ]; then
        # Extract peak names from intersected file (skip header)
        tail -n +2 "${prefix}.crm_intersected.bed" | cut -f4 | sort -u > specified_peaks.txt

        # Create non-intersected file: all peaks EXCEPT those with specified CRM hits
        awk 'BEGIN{OFS="\t"} NR==FNR{specified[\$1]=1; next} !(\$4 in specified)' \\
            specified_peaks.txt $peaks > ${prefix}.crm_non_intersected.bed
    else
        # No specified CRM intersections found - all peaks continue
        cp $peaks ${prefix}.crm_non_intersected.bed
    fi

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


