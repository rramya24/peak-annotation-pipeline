#!/usr/bin/env python3

"""
Extract lncRNA-miRNA relationships from GTF file.
Identifies which lncRNAs encode miRNAs and which miRNAs are encoded by lncRNAs.
"""

import argparse
import sys
import os
import gzip
from collections import defaultdict

def parse_gtf_file(gtf_file):
    """Parse GTF file and extract gene information."""
    genes = {}

    # Handle both gzipped and regular files
    if gtf_file.endswith('.gz'):
        file_handle = gzip.open(gtf_file, 'rt')
    else:
        file_handle = open(gtf_file, 'r')

    try:
        for line in file_handle:
            if line.startswith('#'):
                continue

            line = line.strip()
            if not line:
                continue

            try:
                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                feature_type = parts[2]
                if feature_type != 'gene':
                    continue

                chromosome = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]

                # Parse attributes
                gene_id = None
                gene_name = None
                gene_biotype = None

                # Handle different GTF attribute formats
                attr_pairs = []
                if ';' in attributes:
                    attr_pairs = attributes.split(';')
                else:
                    attr_pairs = [attributes]

                for attr in attr_pairs:
                    attr = attr.strip()
                    if not attr:
                        continue

                    if attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                    elif attr.startswith('gene_name'):
                        gene_name = attr.split('"')[1] if '"' in attr else attr.split()[1]
                    elif attr.startswith('gene_biotype'):
                        gene_biotype = attr.split('"')[1] if '"' in attr else attr.split()[1]
                    elif attr.startswith('gene_type'):
                        gene_biotype = attr.split('"')[1] if '"' in attr else attr.split()[1]

                if gene_id:
                    genes[gene_id] = {
                        'gene_name': gene_name or gene_id,
                        'gene_biotype': gene_biotype or 'unknown',
                        'chromosome': chromosome,
                        'start': start,
                        'end': end,
                        'strand': strand
                    }

            except Exception as e:
                print(f"Warning: Error parsing GTF line: {e}")
                continue

    finally:
        file_handle.close()

    return genes

def find_lncrna_mirna_relationships(genes):
    """Find lncRNA-miRNA encoding relationships."""
    relationships = []

    # Get lncRNAs and miRNAs
    lncrnas = {}
    mirnas = {}

    for gene_id, info in genes.items():
        biotype = info['gene_biotype'].lower()

        if 'lncrna' in biotype or 'lincrna' in biotype or 'long_ncrna' in biotype:
            lncrnas[gene_id] = info
        elif 'mirna' in biotype or 'micro_rna' in biotype or info['gene_name'].startswith('mir-'):
            mirnas[gene_id] = info

    print(f"Found {len(lncrnas)} lncRNAs and {len(mirnas)} miRNAs")

    # Find overlapping relationships
    for lncrna_id, lncrna_info in lncrnas.items():
        for mirna_id, mirna_info in mirnas.items():
            # Check if miRNA is contained within lncRNA
            if (lncrna_info['chromosome'] == mirna_info['chromosome'] and
                lncrna_info['start'] <= mirna_info['start'] and
                lncrna_info['end'] >= mirna_info['end']):

                # Calculate overlap
                overlap_start = max(lncrna_info['start'], mirna_info['start'])
                overlap_end = min(lncrna_info['end'], mirna_info['end'])
                overlap_length = overlap_end - overlap_start + 1
                mirna_length = mirna_info['end'] - mirna_info['start'] + 1
                overlap_fraction = overlap_length / mirna_length

                relationships.append({
                    'lncrna_id': lncrna_id,
                    'lncrna_name': lncrna_info['gene_name'],
                    'mirna_id': mirna_id,
                    'mirna_name': mirna_info['gene_name'],
                    'chromosome': lncrna_info['chromosome'],
                    'overlap_fraction': overlap_fraction,
                    'overlap_info': f"encoded:{overlap_fraction:.2f}"
                })

    print(f"Found {len(relationships)} lncRNA-miRNA encoding relationships")
    return relationships

def write_relationships(relationships, output_file):
    """Write relationships to output file."""
    with open(output_file, 'w') as f:
        f.write("# lncRNA-miRNA encoding relationships\n")
        f.write("# Format: relationship_type\tlncrna_id\tlncrna_name\tlncrna_biotype\tmirna_id\tmirna_name\tmirna_biotype\tchromosome\toverlap_info\n")

        for rel in relationships:
            f.write(f"lncrna_mirna\t{rel['lncrna_id']}\t{rel['lncrna_name']}\tlncRNA\t{rel['mirna_id']}\t{rel['mirna_name']}\tmiRNA\t{rel['chromosome']}\t{rel['overlap_info']}\n")

def write_summary(relationships, genes, summary_file):
    """Write summary statistics."""
    lncrna_count = len([g for g in genes.values() if 'lncrna' in g['gene_biotype'].lower()])
    mirna_count = len([g for g in genes.values() if 'mirna' in g['gene_biotype'].lower() or g['gene_name'].startswith('mir-')])

    # Count unique lncRNAs and miRNAs with relationships
    lncrnas_with_mirnas = set(rel['lncrna_id'] for rel in relationships)
    mirnas_with_lncrnas = set(rel['mirna_id'] for rel in relationships)

    with open(summary_file, 'w') as f:
        f.write("# lncRNA-miRNA Relationship Extraction Summary\n")
        f.write(f"Total genes processed: {len(genes)}\n")
        f.write(f"Total lncRNAs found: {lncrna_count}\n")
        f.write(f"Total miRNAs found: {mirna_count}\n")
        f.write(f"lncRNA-miRNA relationships found: {len(relationships)}\n")
        f.write(f"lncRNAs that encode miRNAs: {len(lncrnas_with_mirnas)}\n")
        f.write(f"miRNAs encoded by lncRNAs: {len(mirnas_with_lncrnas)}\n")
        f.write(f"Average relationships per lncRNA: {len(relationships)/len(lncrnas_with_mirnas):.2f}\n")

        f.write("\n# Relationship details:\n")
        f.write("# lncRNA encodes miRNA - miRNA sequence contained within lncRNA\n")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Extract lncRNA-miRNA encoding relationships from GTF file'
    )
    parser.add_argument(
        '--gtf',
        required=True,
        help='Input GTF file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output relationships file'
    )
    parser.add_argument(
        '--summary',
        required=True,
        help='Output summary file'
    )

    args = parser.parse_args()

    # Check input file
    if not os.path.exists(args.gtf):
        print(f"Error: GTF file not found: {args.gtf}")
        sys.exit(1)

    try:
        print("Parsing GTF file...")
        genes = parse_gtf_file(args.gtf)
        print(f"Parsed {len(genes)} genes from GTF file")

        print("Finding lncRNA-miRNA encoding relationships...")
        relationships = find_lncrna_mirna_relationships(genes)

        print("Writing relationships file...")
        write_relationships(relationships, args.output)

        print("Writing summary...")
        write_summary(relationships, genes, args.summary)

        print("lncRNA-miRNA relationship extraction completed successfully!")

        # Print summary to stdout
        lncrna_count = len([g for g in genes.values() if 'lncrna' in g['gene_biotype'].lower()])
        mirna_count = len([g for g in genes.values() if 'mirna' in g['gene_biotype'].lower() or g['gene_name'].startswith('mir-')])

        print(f"\nSummary:")
        print(f"  Total genes: {len(genes)}")
        print(f"  lncRNAs: {lncrna_count}")
        print(f"  miRNAs: {mirna_count}")
        print(f"  Encoding relationships: {len(relationships)}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

# END OF SCRIPT
