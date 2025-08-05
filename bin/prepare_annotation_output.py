#!/usr/bin/env python3

"""
Prepare final annotation output combining all annotation types with filtering and lncRNA-miRNA expansion.
- Creates cleaned raw output (raw minus problematic HOMER peaks, not genes)
- Filters HOMER peaks that intersect with exons (peak-level filtering)
- Removes targets associated with transposable elements, tRNA, rRNA, and pseudogenes
- Includes lncRNA-miRNA expansion
- Produces multiple output files with different filtering levels
"""

import argparse
import sys
import os
import re
from collections import defaultdict

def parse_gtf_file(gtf_file):
    """Parse GTF file to create gene_id to gene_symbol mapping and extract exon regions."""
    gene_mapping = {}
    gene_biotypes = {}
    exon_regions = set()

    if not os.path.exists(gtf_file):
        print(f"Warning: GTF file not found: {gtf_file}")
        return gene_mapping, gene_biotypes, exon_regions

    print(f"Parsing GTF file: {gtf_file}")

    try:
        with open(gtf_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue

                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                feature_type = parts[2]
                chromosome = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                attributes = parts[8]

                # Extract gene_id
                gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
                if not gene_id_match:
                    continue
                gene_id = gene_id_match.group(1)

                if feature_type == 'gene':
                    # Extract gene_name/gene_symbol
                    gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)
                    gene_symbol = gene_name_match.group(1) if gene_name_match else gene_id

                    # Extract gene_biotype
                    biotype_match = re.search(r'gene_biotype "([^"]+)"', attributes)
                    gene_biotype = biotype_match.group(1) if biotype_match else 'unknown'

                    gene_mapping[gene_id] = gene_symbol
                    gene_biotypes[gene_id] = gene_biotype

                elif feature_type == 'exon':
                    # Store exon regions for filtering
                    exon_regions.add((chromosome, start, end))

                if line_num % 10000 == 0:
                    print(f"  Processed {line_num} lines, found {len(gene_mapping)} genes, {len(exon_regions)} exons")

    except Exception as e:
        print(f"Error parsing GTF file: {e}")
        return gene_mapping, gene_biotypes, exon_regions

    print(f"GTF parsing complete: {len(gene_mapping)} genes, {len(exon_regions)} exon regions found")
    return gene_mapping, gene_biotypes, exon_regions

def parse_homer_file_with_peak_level_filtering(homer_file, exon_regions, gene_mapping, gene_biotypes):
    """Parse HOMER and track problematic peaks (not genes) for removal."""
    genes_raw = {}
    genes_cleaned = {}
    peaks_to_genes = {}
    problematic_peaks = set()
    problematic_genes_count = set()

    if not os.path.exists(homer_file) or os.path.getsize(homer_file) == 0:
        print(f"Warning: HOMER file not found or empty: {homer_file}")
        return genes_raw, genes_cleaned, peaks_to_genes, problematic_peaks

    try:
        with open(homer_file, 'r') as f:
            header = f.readline().strip().split('\t')

            # Find Entrez ID column
            entrez_col = None
            for i, col_name in enumerate(header):
                if 'entrez' in col_name.lower():
                    entrez_col = i
                    break

            if entrez_col is None:
                for pos in [8, 11, 12]:
                    if pos < len(header):
                        entrez_col = pos
                        break

            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 16:
                    peak_id = parts[0]

                    # Extract gene info
                    entrez_id = ''
                    if entrez_col is not None and entrez_col < len(parts):
                        entrez_id = parts[entrez_col].strip()
                        if entrez_id in ['', 'NA', 'N/A', '-', '0']:
                            entrez_id = ''

                    gene_symbol = parts[15] if len(parts) > 15 and parts[15] != '' else 'Unknown'
                    gene_type = parts[17] if len(parts) > 17 and parts[17] != '' else 'unknown'
                    annotation = parts[7] if len(parts) > 7 else 'Unknown'

                    # Determine gene_id
                    if entrez_id:
                        gene_id = entrez_id
                        final_gene_symbol = gene_symbol if gene_symbol != 'Unknown' else entrez_id
                    else:
                        gene_id = gene_symbol
                        final_gene_symbol = gene_symbol

                    if gene_id == 'Unknown' or gene_id == '':
                        continue

                    # Check if this peak should be removed
                    peak_should_remove = False

                    # Check 1: Annotation contains "exon"
                    if 'exon' in annotation.lower():
                        peak_should_remove = True
                        problematic_peaks.add(peak_id)

                    # Check 2: Non-FBgn Entrez ID (like FBtr) - remove entire gene
                    if entrez_id and not entrez_id.startswith('FBgn'):
                        peak_should_remove = True
                        problematic_genes_count.add(gene_id)

                    # Always add to raw collection
                    resolved_biotype = gene_biotypes.get(gene_id, gene_type).lower()
                    if gene_id not in genes_raw:
                        genes_raw[gene_id] = {
                            'gene_id': gene_id,
                            'gene_symbol': final_gene_symbol,
                            'annotation_type': 'homer',
                            'peaks': [],
                            'gene_biotype': resolved_biotype,
                            'homer_annotation': annotation
                        }
                    genes_raw[gene_id]['peaks'].append(peak_id)

                    # Add to cleaned collection only if peak is not problematic AND gene is FBgn
                    if not peak_should_remove or (entrez_id and entrez_id.startswith('FBgn') and 'exon' not in annotation.lower()):
                        if gene_id not in genes_cleaned:
                            genes_cleaned[gene_id] = {
                                'gene_id': gene_id,
                                'gene_symbol': final_gene_symbol,
                                'annotation_type': 'homer',
                                'peaks': [],
                                'gene_biotype': resolved_biotype,
                                'homer_annotation': annotation
                            }

                        # Only add non-problematic peaks
                        if not peak_should_remove:
                            genes_cleaned[gene_id]['peaks'].append(peak_id)

                    peaks_to_genes[peak_id] = gene_id

    except Exception as e:
        print(f"Error parsing HOMER file: {e}")

    # Clean up genes with no remaining peaks in cleaned version
    genes_cleaned_final = {}
    for gene_id, gene_data in genes_cleaned.items():
        if gene_data['peaks']:  # Only keep genes that still have peaks
            genes_cleaned_final[gene_id] = gene_data

    print(f"HOMER: Found {len(genes_raw)} raw genes, {len(problematic_peaks)} problematic peaks removed")
    print(f"HOMER: {len(genes_cleaned_final)} genes remain after peak-level filtering")
    return genes_raw, genes_cleaned_final, peaks_to_genes, problematic_peaks

def parse_homer_file(homer_file, exon_regions=None, filter_exons=False, filter_non_fbgn=False):
    """Parse HOMER annotation file and extract gene information using Entrez ID as gene_id."""
    genes = {}
    peaks_to_genes = {}
    filtered_peaks = 0
    filtered_non_fbgn = 0

    if not os.path.exists(homer_file) or os.path.getsize(homer_file) == 0:
        print(f"Warning: HOMER file not found or empty: {homer_file}")
        return genes, peaks_to_genes

    try:
        with open(homer_file, 'r') as f:
            header = f.readline().strip().split('\t')

            # Find Entrez ID column (common positions: 8, 11, 12)
            entrez_col = None
            for i, col_name in enumerate(header):
                if 'entrez' in col_name.lower():
                    entrez_col = i
                    break

            # If not found by name, try common positions
            if entrez_col is None:
                for pos in [8, 11, 12]:
                    if pos < len(header):
                        entrez_col = pos
                        break

            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 16:  # HOMER format has many columns
                    peak_id = parts[0]

                    # Extract Entrez ID as gene_id
                    entrez_id = ''
                    if entrez_col is not None and entrez_col < len(parts):
                        entrez_id = parts[entrez_col].strip()
                        if entrez_id in ['', 'NA', 'N/A', '-', '0']:
                            entrez_id = ''

                    # Filter non-FBgn IDs if requested
                    if filter_non_fbgn and entrez_id and not entrez_id.startswith('FBgn'):
                        filtered_non_fbgn += 1
                        continue

                    # Extract gene symbol (column 15) and annotation (column 7)
                    gene_symbol = parts[15] if len(parts) > 15 and parts[15] != '' else 'Unknown'
                    gene_type = parts[17] if len(parts) > 17 and parts[17] != '' else 'unknown'
                    annotation = parts[7] if len(parts) > 7 else 'Unknown'

                    # Check if annotation contains "exon" (if filtering enabled)
                    if filter_exons and 'exon' in annotation.lower():
                        filtered_peaks += 1
                        continue  # Skip this peak

                    # Use Entrez ID as gene_id, fallback to gene_symbol if no Entrez ID
                    if entrez_id:
                        gene_id = entrez_id
                        final_gene_symbol = gene_symbol if gene_symbol != 'Unknown' else entrez_id
                    else:
                        gene_id = gene_symbol
                        final_gene_symbol = gene_symbol

                    if gene_id != 'Unknown' and gene_id != '':
                        if gene_id not in genes:
                            genes[gene_id] = {
                                'gene_id': gene_id,
                                'gene_symbol': final_gene_symbol,
                                'annotation_type': 'homer',
                                'peaks': [],
                                'gene_biotype': gene_type,
                                'homer_annotation': annotation
                            }
                        genes[gene_id]['peaks'].append(peak_id)
                        peaks_to_genes[peak_id] = gene_id

    except Exception as e:
        print(f"Error parsing HOMER file: {e}")

    filter_msg = []
    if filter_exons:
        filter_msg.append(f"filtered {filtered_peaks} exon-annotated peaks")
    if filter_non_fbgn:
        filter_msg.append(f"filtered {filtered_non_fbgn} non-FBgn entries")

    filter_status = f" ({', '.join(filter_msg)})" if filter_msg else ""
    print(f"HOMER parsing{filter_status}: Found {len(genes)} genes")

    return genes, peaks_to_genes

def parse_bed_annotation_file(file_path, annotation_type):
    """Parse BED-style annotation files (CRM/intron)."""
    genes = {}
    peaks_to_genes = {}

    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        print(f"Warning: {annotation_type} file not found or empty: {file_path}")
        return genes, peaks_to_genes

    try:
        with open(file_path, 'r') as f:
            # Skip header if present
            header = f.readline().strip()
            if not header.startswith('#') and '\t' in header:
                # Process header line as data
                f.seek(0)

            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) >= 5:
                    gene_id = parts[0] if parts[0] != '' else 'Unknown'
                    gene_symbol = parts[1] if parts[1] != '' else gene_id
                    num_peaks = int(parts[3]) if parts[3].isdigit() else 0
                    peak_names = parts[4] if parts[4] != '' else ''
                    gene_biotype = parts[5] if len(parts) > 5 and parts[5] != '' else 'unknown'

                    if gene_id != 'Unknown' and gene_id != '':
                        if gene_id not in genes:
                            genes[gene_id] = {
                                'gene_id': gene_id,
                                'gene_symbol': gene_symbol,
                                'annotation_type': annotation_type,
                                'peaks': [],
                                'gene_biotype': gene_biotype
                            }

                        if peak_names:
                            peak_list = [p.strip() for p in peak_names.split(',') if p.strip()]
                            genes[gene_id]['peaks'].extend(peak_list)
                            for peak in peak_list:
                                peaks_to_genes[peak] = gene_id

    except Exception as e:
        print(f"Error parsing {annotation_type} file: {e}")

    return genes, peaks_to_genes

def filter_unwanted_biotypes(genes, allowed_biotypes=None):
    """Filter to keep only genes with allowed biotypes (protein_coding and ncRNA)."""
    if allowed_biotypes is None:
        allowed_biotypes = {'protein_coding', 'ncrna', 'nc_rna', 'lncrna', 'mirna', 'snrna', 'snorna'}

    filtered_genes = {}
    removed_count = 0

    for gene_id, gene_data in genes.items():
        biotype = gene_data.get('gene_biotype', '').lower()

        # Check if biotype matches any allowed pattern
        is_allowed = False
        for allowed in allowed_biotypes:
            if allowed in biotype:
                is_allowed = True
                break

        if is_allowed:
            filtered_genes[gene_id] = gene_data
        else:
            removed_count += 1

    print(f"Biotype filtering: Kept only protein_coding and ncRNA genes - Removed {removed_count} genes, kept {len(filtered_genes)} genes")
    return filtered_genes

def parse_lncrna_mirna_mapping(mapping_file):
    """Parse the lncRNA-miRNA mapping file into simple dictionaries."""
    mirna_to_lncrna = {}  # mir-1 -> lncRNA-123
    lncrna_to_mirna = {}  # lncRNA-123 -> [mir-1, mir-2, mir-3]

    if not mapping_file or not os.path.exists(mapping_file):
        print("No lncRNA-miRNA mapping file provided or file not found")
        return mirna_to_lncrna, lncrna_to_mirna

    try:
        with open(mapping_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue

                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    relationship_type = parts[0]
                    gene1_id = parts[1]
                    gene1_name = parts[2]
                    gene2_id = parts[4]
                    gene2_name = parts[5]

                    if relationship_type == 'lncrna_mirna':
                        # lncRNA hosts miRNA
                        mirna_to_lncrna[gene2_id] = {
                            'lncrna_id': gene1_id,
                            'lncrna_name': gene1_name
                        }
                        if gene1_id not in lncrna_to_mirna:
                            lncrna_to_mirna[gene1_id] = []
                        lncrna_to_mirna[gene1_id].append({
                            'mirna_id': gene2_id,
                            'mirna_name': gene2_name
                        })

        print(f"Loaded {len(mirna_to_lncrna)} miRNA->lncRNA mappings and {len(lncrna_to_mirna)} lncRNA->miRNA mappings")
    except Exception as e:
        print(f"Error parsing lncRNA-miRNA mapping file: {e}")

    return mirna_to_lncrna, lncrna_to_mirna

def expand_with_lncrna_mirna(all_gene_ids, gene_mapping, gene_biotypes, mirna_to_lncrna, lncrna_to_mirna):
    """Expand target list with lncRNA-miRNA relationships."""
    if not mirna_to_lncrna and not lncrna_to_mirna:
        print("No lncRNA-miRNA mappings available, skipping expansion")
        return all_gene_ids, []

    expanded_genes = all_gene_ids.copy()
    expansion_log = []
    original_gene_count = len(all_gene_ids)

    # Process each gene in the original target list
    for gene_id, gene_data in list(all_gene_ids.items()):
        # If this target is a miRNA, add its host lncRNA
        if gene_id in mirna_to_lncrna:
            lncrna_info = mirna_to_lncrna[gene_id]
            lncrna_id = lncrna_info['lncrna_id']
            lncrna_name = lncrna_info['lncrna_name']

            if lncrna_id not in expanded_genes:
                # Get proper gene symbol from GTF if available
                proper_symbol = gene_mapping.get(lncrna_id, lncrna_name)
                proper_biotype = gene_biotypes.get(lncrna_id, 'lncRNA')

                expanded_genes[lncrna_id] = {
                    'gene_id': lncrna_id,
                    'gene_symbol': proper_symbol,
                    'annotation_types': [f'lncRNA_encodes_{gene_data["gene_symbol"]}'],
                    'peaks': [f'via_encoded_miRNA_{gene_data["gene_symbol"]}'],
                    'num_peaks': 0,
                    'gene_biotype': proper_biotype
                }
                expansion_log.append(f"Added lncRNA {proper_symbol} ({lncrna_id}) that encodes miRNA {gene_data['gene_symbol']} ({gene_id})")

        # If this target is a lncRNA, add miRNAs it hosts
        if gene_id in lncrna_to_mirna:
            hosted_mirnas = lncrna_to_mirna[gene_id]
            for mirna_info in hosted_mirnas:
                mirna_id = mirna_info['mirna_id']
                mirna_name = mirna_info['mirna_name']

                if mirna_id not in expanded_genes:
                    # Get proper gene symbol from GTF if available
                    proper_symbol = gene_mapping.get(mirna_id, mirna_name)
                    proper_biotype = gene_biotypes.get(mirna_id, 'miRNA')

                    expanded_genes[mirna_id] = {
                        'gene_id': mirna_id,
                        'gene_symbol': proper_symbol,
                        'annotation_types': [f'miRNA_encoded_by_{gene_data["gene_symbol"]}'],
                        'peaks': [f'via_encoding_lncRNA_{gene_data["gene_symbol"]}'],
                        'num_peaks': 0,
                        'gene_biotype': proper_biotype
                    }
                    expansion_log.append(f"Added miRNA {proper_symbol} ({mirna_id}) encoded by lncRNA {gene_data['gene_symbol']} ({gene_id})")

    final_gene_count = len(expanded_genes)
    new_genes_added = final_gene_count - original_gene_count

    print(f"lncRNA-miRNA expansion: {original_gene_count} -> {final_gene_count} genes (+{new_genes_added})")
    return expanded_genes, expansion_log

def combine_annotations(gene_symbol_files, gtf_file, lncrna_mirna_mapping, exon_regions, gene_mapping, gene_biotypes):
    """Combine annotations from all sources with proper merging by gene ID."""

    # Level 1: Raw annotations (no filtering)
    raw_genes = {}

    # Level 2: Cleaned raw (raw with peak-level filtering for HOMER)
    cleaned_raw_genes = {}

    # Level 3: Exon-filtered (HOMER exon filtered + CRM + Intron, merged by gene ID)
    exon_filtered_genes = {}

    all_peaks_to_genes = {}

    for file_path in gene_symbol_files:
        if not file_path or file_path in ["NO_FILE_CRM", "NO_FILE_INTRON", "NO_FILE_HOMER"]:
            continue

        file_name = os.path.basename(file_path).lower()
        print(f"\n=== Processing file: {file_path} ===")

        if 'homer' in file_name or 'annotatepeaks' in file_name:
            # Parse HOMER for raw (no filtering)
            genes_raw, peaks_to_genes_raw = parse_homer_file(file_path, exon_regions, filter_exons=False, filter_non_fbgn=False)

            # Parse HOMER for cleaned raw with peak-level filtering
            genes_raw_for_tracking, genes_cleaned_homer, peaks_to_genes_clean, problematic_peaks = parse_homer_file_with_peak_level_filtering(
                file_path, exon_regions, gene_mapping, gene_biotypes
            )

            # Parse HOMER for exon-filtered (ONLY exon filtering, keep transposons)
            genes_filtered, peaks_to_genes_filtered = parse_homer_file(file_path, exon_regions, filter_exons=True, filter_non_fbgn=False)

            annotation_type = 'homer'
            genes_for_raw = genes_raw
            genes_for_cleaned = genes_cleaned_homer
            genes_for_filtered = genes_filtered

        elif 'crm' in file_name:
            genes, peaks_to_genes = parse_bed_annotation_file(file_path, 'crm')
            annotation_type = 'crm'
            genes_for_raw = genes_for_cleaned = genes_for_filtered = genes
            peaks_to_genes_raw = peaks_to_genes_clean = peaks_to_genes_filtered = peaks_to_genes

        elif 'intron' in file_name:
            genes, peaks_to_genes = parse_bed_annotation_file(file_path, 'intron')
            annotation_type = 'intron'
            genes_for_raw = genes_for_cleaned = genes_for_filtered = genes
            peaks_to_genes_raw = peaks_to_genes_clean = peaks_to_genes_filtered = peaks_to_genes
        else:
            print(f"Unknown file type: {file_path}, skipping...")
            continue

        # Merge into raw_genes (no filtering)
        merge_genes_into_collection(raw_genes, genes_for_raw, annotation_type, gene_mapping, gene_biotypes)

        # Merge into cleaned_raw_genes (peak-level filtering for HOMER)
        merge_genes_into_collection(cleaned_raw_genes, genes_for_cleaned, annotation_type, gene_mapping, gene_biotypes)

        # Merge into exon_filtered_genes (properly combined across all sources)
        merge_genes_into_collection(exon_filtered_genes, genes_for_filtered, annotation_type, gene_mapping, gene_biotypes)

        # Update peak mappings
        all_peaks_to_genes.update(peaks_to_genes_raw if 'peaks_to_genes_raw' in locals() else peaks_to_genes)

    print(f"\nPeak-level filtering applied:")
    print(f"Raw: {len(raw_genes)} genes")
    print(f"Cleaned raw (peak-level filtering): {len(cleaned_raw_genes)} genes")
    print(f"Exon-filtered combined: {len(exon_filtered_genes)} unique genes")

    # Apply lncRNA-miRNA expansion to raw genes
    expanded_raw = {}
    expansion_log_raw = []
    if lncrna_mirna_mapping:
        print("\nPerforming lncRNA-miRNA expansion on raw...")
        mirna_to_lncrna, lncrna_to_mirna = parse_lncrna_mirna_mapping(lncrna_mirna_mapping)
        expanded_raw, expansion_log_raw = expand_with_lncrna_mirna(
            raw_genes, gene_mapping, gene_biotypes, mirna_to_lncrna, lncrna_to_mirna
        )
    else:
        expanded_raw = raw_genes

    # Apply lncRNA-miRNA expansion to cleaned raw
    expanded_cleaned_raw = {}
    expansion_log_cleaned = []
    if lncrna_mirna_mapping:
        print("\nPerforming lncRNA-miRNA expansion on cleaned raw...")
        mirna_to_lncrna, lncrna_to_mirna = parse_lncrna_mirna_mapping(lncrna_mirna_mapping)
        expanded_cleaned_raw, expansion_log_cleaned = expand_with_lncrna_mirna(
            cleaned_raw_genes, gene_mapping, gene_biotypes, mirna_to_lncrna, lncrna_to_mirna
        )
    else:
        expanded_cleaned_raw = cleaned_raw_genes

    # Apply lncRNA-miRNA expansion to exon-filtered genes
    expanded_exon_filtered = {}
    expansion_log_exon = []
    if lncrna_mirna_mapping:
        print("\nPerforming lncRNA-miRNA expansion on exon-filtered...")
        mirna_to_lncrna, lncrna_to_mirna = parse_lncrna_mirna_mapping(lncrna_mirna_mapping)
        expanded_exon_filtered, expansion_log_exon = expand_with_lncrna_mirna(
            exon_filtered_genes, gene_mapping, gene_biotypes, mirna_to_lncrna, lncrna_to_mirna
        )
    else:
        expanded_exon_filtered = exon_filtered_genes

    return raw_genes, expanded_raw, cleaned_raw_genes, expanded_cleaned_raw, exon_filtered_genes, expanded_exon_filtered, all_peaks_to_genes, expansion_log_raw

def merge_genes_into_collection(target_collection, source_genes, annotation_type, gene_mapping, gene_biotypes):
    """Merge genes from source into target collection."""
    for gene_key, gene_data in source_genes.items():
        # Resolve gene symbol from GTF if available
        gene_id = gene_data['gene_id']
        resolved_symbol = gene_mapping.get(gene_id, gene_data['gene_symbol'])
        resolved_biotype = gene_biotypes.get(gene_id, gene_data.get('gene_biotype', 'unknown'))

        if gene_key not in target_collection:
            target_collection[gene_key] = {
                'gene_id': gene_id,
                'gene_symbol': resolved_symbol,
                'annotation_types': [annotation_type],
                'peaks': gene_data['peaks'].copy(),
                'gene_biotype': resolved_biotype
            }
        else:
            # Merge with existing gene
            if annotation_type not in target_collection[gene_key]['annotation_types']:
                target_collection[gene_key]['annotation_types'].append(annotation_type)
            target_collection[gene_key]['peaks'].extend(gene_data['peaks'])

    # Clean up duplicates and add counts
    for gene_key in target_collection:
        target_collection[gene_key]['peaks'] = sorted(list(set(target_collection[gene_key]['peaks'])))
        target_collection[gene_key]['num_peaks'] = len(target_collection[gene_key]['peaks'])
        target_collection[gene_key]['annotation_types'] = sorted(list(set(target_collection[gene_key]['annotation_types'])))

def write_simple_target_file(genes, output_file):
    """Write simple two-column target file."""
    with open(output_file, 'w') as f:
        f.write("gene_id\tgene_symbol\n")
        for gene_key in sorted(genes.keys()):
            gene_data = genes[gene_key]
            f.write(f"{gene_data['gene_id']}\t{gene_data['gene_symbol']}\n")

def write_detailed_target_file(genes, output_file):
    """Write detailed target file."""
    with open(output_file, 'w') as f:
        f.write("gene_id\tgene_symbol\tannotation_types\tnum_peaks\tpeak_names\tgene_biotype\n")
        for gene_key in sorted(genes.keys()):
            gene_data = genes[gene_key]
            annotation_types_str = ','.join(gene_data['annotation_types'])
            peak_names_str = ','.join(gene_data['peaks'])
            f.write(f"{gene_data['gene_id']}\t{gene_data['gene_symbol']}\t{annotation_types_str}\t")
            f.write(f"{gene_data['num_peaks']}\t{peak_names_str}\t{gene_data['gene_biotype']}\n")

def write_summary_file(genes_dict, output_file, prefix):
    """Write summary statistics for different filtering levels."""
    with open(output_file, 'w') as f:
        f.write(f"# Peak Annotation Summary - {prefix}\n")
        f.write("# Multiple filtering levels applied\n\n")

        for level_name, genes in genes_dict.items():
            f.write(f"=== {level_name.upper()} ===\n")
            f.write(f"Total genes: {len(genes)}\n")

            # Count by annotation type
            type_counts = defaultdict(int)
            for gene_data in genes.values():
                for ann_type in gene_data['annotation_types']:
                    type_counts[ann_type] += 1

            f.write("Genes by annotation type:\n")
            for ann_type, count in sorted(type_counts.items()):
                f.write(f"  {ann_type}: {count}\n")

            # Count by biotype
            biotype_counts = defaultdict(int)
            for gene_data in genes.values():
                biotype_counts[gene_data['gene_biotype']] += 1

            f.write("Genes by biotype:\n")
            for biotype, count in sorted(biotype_counts.items()):
                f.write(f"  {biotype}: {count}\n")
            f.write("\n")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Prepare final annotation output with multiple filtering levels and lncRNA-miRNA expansion'
    )
    parser.add_argument(
        '--gene_symbol_files',
        nargs='+',
        required=True,
        help='List of gene symbol annotation files (CRM, intron, HOMER)'
    )
    parser.add_argument(
        '--consensus_peaks',
        help='Consensus peaks file (for reference)'
    )
    parser.add_argument(
        '--gtf_file',
        required=True,
        help='GTF file for gene ID to symbol mapping and exon extraction'
    )
    parser.add_argument(
        '--lncrna_mirna_mapping',
        help='lncRNA-miRNA mapping file (optional)'
    )
    parser.add_argument(
        '--output_prefix',
        required=True,
        help='Output file prefix'
    )

    args = parser.parse_args()

    try:
        print("Starting final annotation output preparation with peak-level filtering...")
        print(f"Input files: {args.gene_symbol_files}")
        print(f"GTF file: {args.gtf_file}")
        print(f"Output prefix: {args.output_prefix}")

        # Parse GTF file
        gene_mapping, gene_biotypes, exon_regions = parse_gtf_file(args.gtf_file)

        # Combine all annotations with different filtering levels
        raw_genes, expanded_raw, cleaned_raw_genes, expanded_cleaned_raw, exon_filtered_genes, expanded_exon_filtered, all_peaks_to_genes, expansion_log_raw = combine_annotations(
            args.gene_symbol_files, args.gtf_file, args.lncrna_mirna_mapping,
            exon_regions, gene_mapping, gene_biotypes
        )

        # Write output files for each filtering level (matching Nextflow expectations)

        # 1. Raw targets (no filtering)
        write_simple_target_file(raw_genes, f"{args.output_prefix}.raw_targets.tsv")
        write_detailed_target_file(raw_genes, f"{args.output_prefix}.raw_targets_detailed.tsv")

        # 2. Raw expanded targets (raw + lncRNA expansion)
        write_simple_target_file(expanded_raw, f"{args.output_prefix}.raw_expanded_targets.tsv")
        write_detailed_target_file(expanded_raw, f"{args.output_prefix}.raw_expanded_targets_detailed.tsv")

        # 3. Cleaned raw targets (peak-level filtering for HOMER)
        write_simple_target_file(cleaned_raw_genes, f"{args.output_prefix}.cleaned_raw_targets.tsv")
        write_detailed_target_file(cleaned_raw_genes, f"{args.output_prefix}.cleaned_raw_targets_detailed.tsv")

        # 4. Final cleaned + expanded targets (cleaned raw + lncRNA expansion) - MAIN OUTPUT
        write_simple_target_file(expanded_cleaned_raw, f"{args.output_prefix}.final_cleaned_targets.tsv")
        write_detailed_target_file(expanded_cleaned_raw, f"{args.output_prefix}.final_cleaned_targets_detailed.tsv")

        # 5. Exon-filtered targets
        write_simple_target_file(exon_filtered_genes, f"{args.output_prefix}.exon_filtered_targets.tsv")
        write_detailed_target_file(exon_filtered_genes, f"{args.output_prefix}.exon_filtered_targets_detailed.tsv")

        # 6. Final putative targets (exon-filtered + lncRNA expansion)
        write_simple_target_file(expanded_exon_filtered, f"{args.output_prefix}.final_putative_targets.tsv")
        write_detailed_target_file(expanded_exon_filtered, f"{args.output_prefix}.final_putative_targets_detailed.tsv")

        # 7. Legacy all_target_genes file (for compatibility) - using final cleaned targets
        write_detailed_target_file(expanded_cleaned_raw, f"{args.output_prefix}.all_target_genes.txt")

        # 8. Summary file
        genes_dict = {
            'raw': raw_genes,
            'raw_expanded': expanded_raw,
            'cleaned_raw': cleaned_raw_genes,
            'final_cleaned_expanded': expanded_cleaned_raw,
            'exon_filtered': exon_filtered_genes,
            'final_putative_expanded': expanded_exon_filtered
        }
        write_summary_file(genes_dict, f"{args.output_prefix}.filtering_summary.txt", args.output_prefix)

        # 9. Expansion log for raw
        if expansion_log_raw:
            with open(f"{args.output_prefix}.lncrna_mirna_expansion_raw.log", 'w') as f:
                f.write("# lncRNA-miRNA Target Expansion Log (Raw)\n")
                f.write(f"# Total expansions: {len(expansion_log_raw)}\n")
                f.write("# Expansion details:\n")
                for log_entry in expansion_log_raw:
                    f.write(f"{log_entry}\n")

        print(f"\nFinal annotation output completed successfully!")
        print(f"Summary:")
        print(f"  Raw targets: {len(raw_genes)}")
        print(f"  Raw expanded targets: {len(expanded_raw)} (raw + lncRNA expansion)")
        print(f"  Cleaned raw targets: {len(cleaned_raw_genes)} (peak-level filtering)")
        print(f"  Final cleaned + expanded: {len(expanded_cleaned_raw)} (your main output)")
        print(f"  Exon-filtered targets: {len(exon_filtered_genes)}")
        print(f"  Final putative targets: {len(expanded_exon_filtered)} (exon-filtered + lncRNA expansion)")
        if expansion_log_raw:
            print(f"  lncRNA-miRNA expansions: {len(expansion_log_raw)}")
        print(f"  Output files created with prefix: {args.output_prefix}")
        print(f"\n** Your main output file: {args.output_prefix}.final_cleaned_targets.tsv **")
        print(f"** Gene FBgn0050009 should now be retained in final_cleaned_targets.tsv **")

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
