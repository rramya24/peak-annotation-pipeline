# Test Data for Peak Annotation Pipeline

This directory contains test data files for the multi-step peak annotation pipeline.

## Files

### Peak Files (peaks/)
- `rep1_peaks.bed` - Replicate 1 peak calls
- `rep2_peaks.bed` - Replicate 2 peak calls
- `rep3_peaks.bed` - Replicate 3 peak calls

### Annotation Files (annotations/)
- `crm_regions.bed` - Cis-regulatory module coordinates no chr and gene-names converted to fbgn geneids
- `first_intron_regions.bed` - First intron coordinates no chr and gene name converted to geneids
- `drosophila_ensembl_114.gtf` - drosophila gtf file ensembl version 114
- `Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa` - the corresponding fasta file dm6 Ensembl version 114

## Usage

These files are used by the test configurations:
- `conf/test.config` - Quick test with limited resources
- `conf/test_full.config` - Full test with all features


