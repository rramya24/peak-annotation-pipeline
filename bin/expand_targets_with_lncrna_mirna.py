#!/usr/bin/env python3

"""
Simple lncRNA-miRNA expansion for putative targets:
- If peak hits miRNA → add its host lncRNA to final target list
- If peak hits lncRNA → add miRNAs it encodes to final target list
"""

import argparse
import os
import sys

def parse_lncrna_mirna_mapping(mapping_file):
    """Parse the lncRNA-miRNA mapping file into simple dictionaries."""
    mirna_to_lncrna = {}  # mir-1 -> lncRNA-123
    lncrna_to_mirna = {}  # lncRNA-123 -> [mir-1, mir-2, mir-3]

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

    return mirna_to_lncrna, lncrna_to_mirna

def parse_target_genes(target_file):
    """Parse target genes file."""
    targets = []

    with open(target_file, 'r') as f:
        header = f.readline().strip()
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) >= 2:
                gene_id = parts[0]
                gene_symbol = parts[1]
                annotation_types = parts[2] if len(parts) > 2 else ''
                num_peaks = parts[3] if len(parts) > 3 else '0'
                peak_names = parts[4] if len(parts) > 4 else ''

                targets.append({
                    'gene_id': gene_id,
                    'gene_symbol': gene_symbol,
                    'annotation_types': annotation_types,
                    'num_peaks': num_peaks,
                    'peak_names': peak_names,
                    'is_original': True
                })

    return targets

def simple_expansion(targets, mirna_to_lncrna, lncrna_to_mirna):
    """Simple expansion: add direct relationships only."""
    expanded_targets = []
    added_genes = set()
    expansion_log = []

    # First, add all original targets
    for target in targets:
        expanded_targets.append(target)
        added_genes.add(target['gene_id'])

    # Then, add simple relationships
    for target in targets:
        gene_id = target['gene_id']
        gene_symbol = target['gene_symbol']

        # If this target is a miRNA, add its host lncRNA
        if gene_id in mirna_to_lncrna:
            lncrna_info = mirna_to_lncrna[gene_id]
            lncrna_id = lncrna_info['lncrna_id']
            lncrna_name = lncrna_info['lncrna_name']

            if lncrna_id not in added_genes:
                expanded_targets.append({
                    'gene_id': lncrna_id,
                    'gene_symbol': lncrna_name,
                    'annotation_types': f'lncRNA_encodes_{gene_symbol}',
                    'num_peaks': '0',
                    'peak_names': f'via_encoded_miRNA_{gene_symbol}',
                    'is_original': False
                })
                added_genes.add(lncrna_id)
                expansion_log.append(f"Added lncRNA {lncrna_name} that encodes miRNA {gene_symbol}")

        # If this target is a lncRNA, add miRNAs it hosts
        if gene_id in lncrna_to_mirna:
            hosted_mirnas = lncrna_to_mirna[gene_id]
            for mirna_info in hosted_mirnas:
                mirna_id = mirna_info['mirna_id']
                mirna_name = mirna_info['mirna_name']

                if mirna_id not in added_genes:
                    expanded_targets.append({
                        'gene_id': mirna_id,
                        'gene_symbol': mirna_name,
                        'annotation_types': f'miRNA_encoded_by_{gene_symbol}',
                        'num_peaks': '0',
                        'peak_names': f'via_encoding_lncRNA_{gene_symbol}',
                        'is_original': False
                    })
                    added_genes.add(mirna_id)
                    expansion_log.append(f"Added miRNA {mirna_name} encoded by lncRNA {gene_symbol}")

    return expanded_targets, expansion_log

def write_expanded_targets(expanded_targets, output_file):
    """Write expanded targets to output file."""
    with open(output_file, 'w') as f:
        # Write header
        f.write("gene_id\tgene_symbol\tannotation_types\tnum_peaks\tpeak_names\tis_original\n")

        for target in expanded_targets:
            gene_id = target['gene_id']
            gene_symbol = target['gene_symbol']
            annotation_types = target['annotation_types']
            num_peaks = target['num_peaks']
            peak_names = target['peak_names']
            is_original = 'TRUE' if target['is_original'] else 'FALSE'

            f.write(f"{gene_id}\t{gene_symbol}\t{annotation_types}\t{num_peaks}\t{peak_names}\t{is_original}\n")

def write_expansion_log(expansion_log, log_file):
    """Write expansion log to file."""
    with open(log_file, 'w') as f:
        f.write("# Simple lncRNA-miRNA Target Expansion Log\n")
        f.write(f"# Total expansions: {len(expansion_log)}\n")
        f.write("# Expansion details:\n")

        for log_entry in expansion_log:
            f.write(f"{log_entry}\n")

def generate_expansion_summary(original_targets, expanded_targets, summary_file):
    """Generate summary of expansion results."""
    original_count = len(original_targets)
    expanded_count = len(expanded_targets)
    new_count = expanded_count - original_count

    # Count additions by type
    lncrna_hosts = 0
    hosted_mirnas = 0

    for target in expanded_targets:
        if not target['is_original']:
            if 'lncRNA_encodes' in target['annotation_types']:
                lncrna_hosts += 1
            elif 'miRNA_encoded_by' in target['annotation_types']:
                hosted_mirnas += 1

    with open(summary_file, 'w') as f:
        f.write("# Simple lncRNA-miRNA Target Expansion Summary\n")
        f.write(f"Original targets: {original_count}\n")
        f.write(f"Expanded targets: {expanded_count}\n")
        f.write(f"New targets added: {new_count}\n")

        # Fix potential division by zero
        if original_count > 0:
            expansion_ratio = expanded_count / original_count
            f.write(f"Expansion ratio: {expansion_ratio:.2f}\n")
        else:
            f.write("Expansion ratio: 0.00 (no original targets)\n")

        f.write("\n# New targets by type:\n")
        f.write(f"lncRNA hosts added: {lncrna_hosts}\n")
        f.write(f"Hosted miRNAs added: {hosted_mirnas}\n")
        f.write("\n# Simple expansion logic:\n")
        f.write("- If peak hits miRNA → add its host lncRNA\n")
        f.write("- If peak hits lncRNA → add miRNAs it hosts\n")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Simple lncRNA-miRNA expansion for putative targets'
    )
    parser.add_argument(
        '--target-genes',
        required=True,
        help='Input target genes file'
    )
    parser.add_argument(
        '--lncrna-mirna-mapping',
        required=True,
        help='lncRNA-miRNA mapping file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output expanded targets file'
    )
    parser.add_argument(
        '--log',
        required=True,
        help='Output expansion log file'
    )
    parser.add_argument(
        '--summary',
        required=True,
        help='Output expansion summary file'
    )

    args = parser.parse_args()

    # Validate input files
    if not os.path.exists(args.target_genes):
        print(f"Error: Target genes file not found: {args.target_genes}")
        sys.exit(1)

    if not os.path.exists(args.lncrna_mirna_mapping):
        print(f"Error: lncRNA-miRNA mapping file not found: {args.lncrna_mirna_mapping}")
        sys.exit(1)

    try:
        # Parse input files
        print("Parsing lncRNA-miRNA relationships...")
        mirna_to_lncrna, lncrna_to_mirna = parse_lncrna_mirna_mapping(args.lncrna_mirna_mapping)

        print("Parsing target genes...")
        targets = parse_target_genes(args.target_genes)

        print(f"Found {len(targets)} original targets")

        # Simple expansion
        print("Performing simple lncRNA-miRNA expansion...")
        expanded_targets, expansion_log = simple_expansion(targets, mirna_to_lncrna, lncrna_to_mirna)

        print(f"Expanded to {len(expanded_targets)} targets ({len(expansion_log)} additions)")

        # Write output files
        print("Writing expanded targets...")
        write_expanded_targets(expanded_targets, args.output)

        print("Writing expansion log...")
        write_expansion_log(expansion_log, args.log)

        print("Writing expansion summary...")
        generate_expansion_summary(targets, expanded_targets, args.summary)

        print("Simple lncRNA-miRNA expansion completed successfully!")

        # Print summary to stdout
        original_count = len(targets)
        expanded_count = len(expanded_targets)
        new_count = expanded_count - original_count

        print(f"\nSummary:")
        print(f"  Original targets: {original_count}")
        print(f"  Expanded targets: {expanded_count}")
        print(f"  New targets added: {new_count}")
        print(f"  Expansion ratio: {expanded_count/original_count:.2f}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

# END OF SCRIPT
