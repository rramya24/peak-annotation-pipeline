class GenomeSpeciesMapping {

    // Map genome names to Ensembl species names
    static def genomeToSpecies = [
        'dm6': 'drosophila_melanogaster',
        'dm3': 'drosophila_melanogaster',
        'hg38': 'homo_sapiens',
        'hg19': 'homo_sapiens',
        'GRCh38': 'homo_sapiens',
        'GRCh37': 'homo_sapiens',
        'mm10': 'mus_musculus',
        'mm9': 'mus_musculus',
        'GRCm38': 'mus_musculus',
        'GRCm37': 'mus_musculus',
        'rn6': 'rattus_norvegicus',
        'rn5': 'rattus_norvegicus',
        'ce11': 'caenorhabditis_elegans',
        'ce10': 'caenorhabditis_elegans',
        'danRer11': 'danio_rerio',
        'danRer10': 'danio_rerio',
        'galGal6': 'gallus_gallus',
        'galGal5': 'gallus_gallus',
        'susScr11': 'sus_scrofa',
        'susScr3': 'sus_scrofa',
        'bosTau9': 'bos_taurus',
        'bosTau8': 'bos_taurus',
        'canFam3': 'canis_familiaris',
        'canFam6': 'canis_familiaris'
    ]

    // Get species name from genome
    static def getSpeciesFromGenome(genome) {
        return genomeToSpecies.get(genome?.toLowerCase(), genome)
    }

    // Get GTF URL for species and version
    static def getGtfUrl(species, version) {
        return "https://ftp.ensembl.org/pub/release-${version}/gtf/${species}/"
    }

    // Get common genome information
    static def getGenomeInfo(genome) {
        def species = getSpeciesFromGenome(genome)
        def info = [
            genome: genome,
            species: species,
            ensembl_name: species,
            homer_name: genome
        ]

        // Add specific information for common genomes
        switch(genome?.toLowerCase()) {
            case 'dm6':
                info.common_name = 'Drosophila melanogaster'
                info.assembly = 'BDGP6'
                break
            case 'hg38':
                info.common_name = 'Homo sapiens'
                info.assembly = 'GRCh38'
                break
            case 'mm10':
                info.common_name = 'Mus musculus'
                info.assembly = 'GRCm38'
                break
            case 'rn6':
                info.common_name = 'Rattus norvegicus'
                info.assembly = 'Rnor_6.0'
                break
            case 'ce11':
                info.common_name = 'Caenorhabditis elegans'
                info.assembly = 'WBcel235'
                break
            case 'danrer11':
                info.common_name = 'Danio rerio'
                info.assembly = 'GRCz11'
                break
            default:
                info.common_name = species?.replaceAll('_', ' ')?.split(' ')?.collect { it.capitalize() }?.join(' ')
                info.assembly = 'Unknown'
        }

        return info
    }

    // Check if genome is supported
    static def isGenomeSupported(genome) {
        return genomeToSpecies.containsKey(genome?.toLowerCase())
    }

    // Get list of supported genomes
    static def getSupportedGenomes() {
        return genomeToSpecies.keySet().sort()
    }

    // Get help text for genome parameter
    static def getGenomeHelpText() {
        def supported = getSupportedGenomes().join(', ')
        return """
        Supported genomes: ${supported}

        For other genomes, please provide:
        - --species: Ensembl species name (e.g., 'drosophila_melanogaster')
        - --gtf: Direct path to GTF file

        Example custom usage:
        --genome custom_genome --species your_species_name --gtf /path/to/your.gtf
        """
    }
}

