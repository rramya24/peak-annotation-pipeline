process DOWNLOAD_GTF {
    tag "${species}_${ensembl_version}"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::wget=1.20.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    val species
    val ensembl_version

    output:
    path "*.gtf.gz", emit: gtf
    path "*.download_info.txt", emit: info
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Download GTF from Ensembl
    echo "Downloading GTF for ${species} from Ensembl release ${ensembl_version}"

    # Construct Ensembl FTP URL
    ENSEMBL_URL="https://ftp.ensembl.org/pub/release-${ensembl_version}/gtf/${species}/"

    # Get the GTF filename
    GTF_FILE=\$(curl -s \$ENSEMBL_URL | grep -o "${species}.*\\.gtf\\.gz" | head -1)

    if [ -z "\$GTF_FILE" ]; then
        echo "Error: Could not find GTF file for ${species} in Ensembl release ${ensembl_version}"
        echo "Available species can be found at: https://ftp.ensembl.org/pub/release-${ensembl_version}/gtf/"
        exit 1
    fi

    # Download the GTF file
    echo "Downloading: \$GTF_FILE"
    wget \${ENSEMBL_URL}\$GTF_FILE

    # Create download info file
    echo "Species: ${species}" > \${GTF_FILE%.gtf.gz}.download_info.txt
    echo "Ensembl version: ${ensembl_version}" >> \${GTF_FILE%.gtf.gz}.download_info.txt
    echo "Download URL: \${ENSEMBL_URL}\$GTF_FILE" >> \${GTF_FILE%.gtf.gz}.download_info.txt
    echo "Download date: \$(date)" >> \${GTF_FILE%.gtf.gz}.download_info.txt
    echo "File size: \$(ls -lh \$GTF_FILE | awk '{print \$5}')" >> \${GTF_FILE%.gtf.gz}.download_info.txt

    # Verify download
    if [ ! -f "\$GTF_FILE" ]; then
        echo "Error: Download failed"
        exit 1
    fi

    echo "GTF download completed: \$GTF_FILE"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | sed 's/GNU Wget //g')
        ensembl: ${ensembl_version}
    END_VERSIONS
    """

    stub:
    """
    touch ${species}.${ensembl_version}.gtf.gz
    touch ${species}.${ensembl_version}.download_info.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | sed 's/GNU Wget //g')
        ensembl: ${ensembl_version}
    END_VERSIONS
    """
}



