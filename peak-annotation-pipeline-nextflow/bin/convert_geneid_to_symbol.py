#!/usr/bin/env python3

"""
Convert gene IDs to gene symbols using GTF annotation.
Preserves peak names throughout the process.
"""

import argparse
import sys
import os
import gzip

def parse_gtf_file(gtf_file):
    """Parse GTF file and extract gene ID to symbol mapping."""
    gene_mapping = {}

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
                    gene_mapping[gene_id] = {
                        'gene_symbol': gene_name or gene_id,
                        'gene_biotype': gene_biotype or 'unknown'
                    }

            except Exception as e:
                print(f"Warning: Error parsing GTF line: {e}")
                continue

    finally:
        file_handle.close()

    return gene_mapping

def convert_gene_ids(annotation_file, gtf_file, output_file):
    """Convert gene IDs to symbols and preserve peak names."""
    try:
        # Parse GTF file
        print("Parsing GTF file for gene mappings...")
        gene_mapping = parse_gtf_file(gtf_file)
        print(f"Found {len(gene_mapping)} gene mappings")

        # Read annotation file
        annotations = []
        with open(annotation_file, 'r') as f:
            header = f.readline().strip()
            header_parts = header.split('\t')

            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 4:
                    gene_id = parts[0]
                    old_gene_symbol = parts[1]
                    annotation_type = parts[2]
                    num_peaks = parts[3]
                    peak_names = parts[4] if len(parts) > 4 else ''

                    # Additional columns (for CRM/intron regions)
                    extra_columns = parts[5:] if len(parts) > 5 else []

                    # Convert gene ID to symbol
                    if gene_id in gene_mapping:
                        new_gene_symbol = gene_mapping[gene_id]['gene_symbol']
                        gene_biotype = gene_mapping[gene_id]['gene_biotype']
                    else:
                        # Try to extract from existing symbol or use as-is
                        new_gene_symbol = old_gene_symbol
                        gene_biotype = 'unknown'

                        # Handle special cases for CRM and intron annotations
                        if gene_id.startswith('CRM_'):
                            new_gene_symbol = gene_id.replace('CRM_', '')
                            gene_biotype = 'crm_target'
                        elif gene_id.startswith('INTRON_'):
                            new_gene_symbol = gene_id.replace('INTRON_', '')
                            gene_biotype = 'intron_target'

                    annotations.append({
                        'gene_id': gene_id,
                        'gene_symbol': new_gene_symbol,
                        'annotation_type': annotation_type,
                        'num_peaks': num_peaks,
                        'peak_names': peak_names,
                        'gene_biotype': gene_biotype,
                        'extra_columns': extra_columns
                    })

        # Write converted annotations
        with open(output_file, 'w') as f:
            # Write header - preserve original structure but add gene_biotype
            if 'gene_biotype' not in header:
                header_parts_new = header_parts[:5] + ['gene_biotype'] + header_parts[5:]
                f.write('\t'.join(header_parts_new) + '\n')
            else:
                f.write(header + '\n')

            for annotation in annotations:
                # Write main columns
                f.write(f"{annotation['gene_id']}\t{annotation['gene_symbol']}\t{annotation['annotation_type']}\t{annotation['num_peaks']}\t{annotation['peak_names']}\t{annotation['gene_biotype']}")

                # Write extra columns if present
                if annotation['extra_columns']:
                    f.write('\t' + '\t'.join(annotation['extra_columns']))

                f.write('\n')

        print(f"Converted {len(annotations)} gene annotations")
        print(f"Results written to: {output_file}")

        # Print summary
        converted_count = sum(1 for a in annotations if a['gene_id'] in gene_mapping)
        print(f"Successfully converted: {converted_count}/{len(annotations)} genes")

        return True

    except Exception as e:
        print(f"Error in gene ID conversion: {e}")
        return False

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Convert gene IDs to gene symbols using GTF annotation'
    )
    parser.add_argument(
        '--annotation',
        required=True,
        help='Input annotation file with gene IDs'
    )
    parser.add_argument(
        '--gtf',
        required=True,
        help='GTF annotation file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output file with gene symbols'
    )

    args = parser.parse_args()

    # Check input files
    if not os.path.exists(args.annotation):
        print(f"Error: Annotation file not found: {args.annotation}")
        sys.exit(1)

    if not os.path.exists(args.gtf):
        print(f"Error: GTF file not found: {args.gtf}")
        sys.exit(1)

    try:
        print("Converting gene IDs to symbols...")
        success = convert_gene_ids(args.annotation, args.gtf, args.output)

        if not success:
            print("Error: Failed to convert gene IDs")
            sys.exit(1)

        print("Gene ID to symbol conversion completed successfully!")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

# END OF SCRIPT
