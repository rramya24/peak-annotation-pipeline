# File Implementation Checklist

## Core Pipeline Files
- [ ] main.nf
- [ ] nextflow.config
- [ ] workflows/multistep_peak_annotation.nf

## Python Scripts (bin/)
- [ ] check_samplesheet.py
- [ ] consensus_peaks.py
- [ ] intersect_crm_annotate.py
- [ ] intersect_intron_annotate.py
- [ ] convert_geneid_to_symbol.py
- [ ] prepare_annotation_output.py
- [ ] extract_first_introns.py
- [ ] extract_lncrna_mirna_relationships.py
- [ ] expand_targets_with_lncrna_mirna.py

## Configuration Files (conf/)
- [ ] base.config
- [ ] modules.config
- [ ] test.config

## Library Files (lib/)
- [ ] GenomeSpeciesMapping.groovy

## Local Modules (modules/local/)
- [ ] bedtools_intersect_crm/main.nf
- [ ] bedtools_intersect_intron/main.nf
- [ ] consensus_peaks/main.nf
- [ ] convert_geneid_to_symbol/main.nf
- [ ] prepare_annotation_output/main.nf
- [ ] samplesheet_check/main.nf
- [ ] download_gtf/main.nf
- [ ] extract_introns/main.nf
- [ ] extract_lncrna_mirna/main.nf
- [ ] expand_targets_lncrna/main.nf

## nf-core Modules (modules/nf-core/)
- [ ] gunzip/main.nf
- [ ] homer/annotatepeaks/main.nf

## Subworkflows (subworkflows/local/)
- [ ] input_check.nf
- [ ] macs2_consensus.nf
- [ ] peak_annotation.nf

## Documentation
- [ ] README.md

## Standard nf-core Files (use templates)
- [ ] .nf-core.yml
- [ ] .gitignore
- [ ] CHANGELOG.md
- [ ] CITATIONS.md
- [ ] CODE_OF_CONDUCT.md
- [ ] LICENSE
- [ ] nextflow_schema.json
- [ ] modules/nf-core/custom/dumpsoftwareversions/main.nf
- [ ] modules/nf-core/multiqc/main.nf
- [ ] lib/NfcoreTemplate.groovy
- [ ] lib/WorkflowMain.groovy

## Testing
- [ ] Test with minimal dataset
- [ ] Test on UCL Myriad
- [ ] Validate all outputs

