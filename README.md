# multi-step peak annotation pipeline
This uses netxflow workflow manager to perform a multi-step peak annotation for peaks obtained from ChIP seq or DamID seq data.
The pipeline currently is made to annotated Drosophila melanogaster data at dm6 (based on the test sample data given) but can be adapted to any organism if given the correct input files
The peaks in bed format is first annotated to known Cis Regulatory Motifs (CRMs) for Droosphila the data was pulled out from redfly database using intersect bed 
The peaks are assigned to genes if it is found within its 1st intron (as most regulatory regions are within the 1st introns of genes
The peaks are also annotated using homer: assigning peaks to the genes with the most proximal TSS
The results has a detail annotation of where the peaks are located so that the annotation of the peak to a gene can be in the order preference currently set at CRM >1st intron > TSS
