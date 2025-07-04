#!/usr/bin/env python3

"""
Extract first introns from GTF file.
Creates BED file with first intron coordinates for each gene.
"""

import argparse
import os
import sys
import re
from collections import defaultdict

def parse_gtf_attributes(attribute_string):
    """Parse GTF attribute string into dictionary."""
    attributes = {}
    # Split by semicolon and process each attribute
    for attr in attribute_string.split(';'):
        attr = attr.strip()
        if attr:
            # Match attribute name and value
            match = re.match(r'(\w+)\s*["\s]*([^"]+)["\s]*', attr)
            if match:
                key, value = match.groups()
                attributes[key] = value.strip('"')
    return attributes

def extract_gene_structure(gtf_file):
    """Extract gene structure from GTF file."""
    genes = defaultdict(lambda: {
        'chrom': '',
        'start': float('inf'),
        'end': 0,
        'strand': '',
        'gene_name': '',
        'exons': []
    })

    print(f"Reading GTF file: {gtf_file}")

    with open(gtf_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue

            line = line.strip()
            if not line:
                continue

            try:
                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                chrom = parts[0]
                feature = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = parse_gtf_attributes(parts[8])

                # Process genes and exons
                if feature == 'gene':
                    gene_id = attributes.get('gene_id', '')
                    gene_name = attributes.get('gene_name', gene_id)

                    if gene_id:
                        genes[gene_id]['chrom'] = chrom
                        genes[gene_id]['start'] = start
                        genes[gene_id]['end'] = end
                        genes[gene_id]['strand'] = strand
                        genes[gene_id]['gene_name'] = gene_name

                elif feature == 'exon':
                    gene_id = attributes.get('gene_id', '')
                    if gene_id:
                        genes[gene_id]['exons'].append({
                            'start': start,
                            'end': end,
                            'transcript_id': attributes.get('transcript_id', ''),
                            'exon_number': attributes.get('exon_number', '1')
                        })

            except Exception as e:
                print(f"Warning: Error parsing line {line_num}: {e}")
                continue

    print(f"Processed {len(genes)} genes")
    return genes

def find_first_introns(genes):
    """Find first introns for each gene."""
    first_introns = []

    for gene_id, gene_info in genes.items():
        if len(gene_info['exons']) < 2:
            continue  # Need at least 2 exons to have an intron

        # Group exons by transcript
        transcript_exons = defaultdict(list)
        for exon in gene_info['exons']:
            transcript_id = exon['transcript_id']
            if transcript_id:
                transcript_exons[transcript_id].append(exon)

        # Find first intron for each transcript
        transcript_first_introns = []

        for transcript_id, exons in transcript_exons.items():
            if len(exons) < 2:
                continue

            # Sort exons by position
            exons.sort(key=lambda x: x['start'])

            # Find first intron (between first and second exon)
            if gene_info['strand'] == '+':
                # Forward strand: first intron is between exon 1 and exon 2
                intron_start = exons[0]['end'] + 1
                intron_end = exons[1]['start'] - 1
            else:
                # Reverse strand: first intron is between last and second-to-last exon
                intron_start = exons[-2]['end'] + 1
                intron_end = exons[-1]['start'] - 1

            # Ensure valid coordinates
            if intron_start <= intron_end:
                transcript_first_introns.append({
                    'start': intron_start,
                    'end': intron_end,
                    'transcript_id': transcript_id
                })

        # Use the longest first intron as representative
        if transcript_first_introns:
            # Sort by length (longest first)
            transcript_first_introns.sort(key=lambda x: x['end'] - x['start'], reverse=True)
            best_intron = transcript_first_introns[0]

            first_introns.append({
                'chrom': gene_info['chrom'],
                'start': best_intron['start'],
                'end': best_intron['end'],
                'gene_id': gene_id,
                'gene_name': gene_info['gene_name'],
                'strand': gene_info['strand'],
                'transcript_id': best_intron['transcript_id']
            })

    return first_introns

def write_introns_bed(first_introns, output_file):
    """Write first introns to BED file."""
    print(f"Writing {len(first_introns)} first introns to {output_file}")

    with open(output_file, 'w') as f:
        f.write("# First introns extracted from GTF\n")
        f.write("# chrom\tstart\tend\tgene_id\tscore\tstrand\n")

        for intron in sorted(first_introns, key=lambda x: (x['chrom'], x['start'])):
            f.write(f"{intron['chrom']}\t{intron['start']}\t{intron['end']}\t{intron['gene_id']}\t0\t{intron['strand']}\n")

def write_log(genes, first_introns, gtf_file, output_file, log_file):
    """Write summary log."""
    with open(log_file, 'w') as f:
        f.write("First Intron Extraction Summary\n")
        f.write("="*40 + "\n")
        f.write(f"Input GTF: {gtf_file}\n")
        f.write(f"Output BED: {output_file}\n")
        f.write(f"Total genes processed: {len(genes)}\n")
        f.write(f"Genes with first introns: {len(first_introns)}\n")
        f.write(f"Genes without introns: {len(genes) - len(first_introns)}\n")
        f.write("\n")

        # Count by chromosome
        chrom_counts = defaultdict(int)
        for intron in first_introns:
            chrom_counts[intron['chrom']] += 1

        f.write("First introns by chromosome:\n")
        for chrom in sorted(chrom_counts.keys()):
            f.write(f"  {chrom}: {chrom_counts[chrom]} introns\n")

        f.write("\n")

        # Length statistics
        lengths = [intron['end'] - intron['start'] + 1 for intron in first_introns]
        if lengths:
            f.write(f"Intron length statistics:\n")
            f.write(f"  Min length: {min(lengths)} bp\n")
            f.write(f"  Max length: {max(lengths)} bp\n")
            f.write(f"  Mean length: {sum(lengths)/len(lengths):.1f} bp\n")

def main():
    parser = argparse.ArgumentParser(description='Extract first introns from GTF file')
    parser.add_argument('--gtf', required=True, help='Input GTF file')
    parser.add_argument('--output', required=True, help='Output BED file')
    parser.add_argument('--log', required=True, help='Log file')

    args = parser.parse_args()

    # Check input file
    if not os.path.exists(args.gtf):
        print(f"Error: GTF file {args.gtf} does not exist")
        sys.exit(1)

    # Extract gene structure
    print("Extracting gene structure from GTF...")
    genes = extract_gene_structure(args.gtf)

    # Find first introns
    print("Finding first introns...")
    first_introns = find_first_introns(genes)

    # Write output
    write_introns_bed(first_introns, args.output)
    write_log(genes, first_introns, args.gtf, args.output, args.log)

    print(f"First intron extraction complete:")
    print(f"  Total genes: {len(genes)}")
    print(f"  Genes with first introns: {len(first_introns)}")
    print(f"  Output: {args.output}")

if __name__ == "__main__":
    main()

# END OF SCRIPT
