#!/usr/bin/env python3

"""
Intersect peaks with CRM regions and annotate with gene information.
FILTERS OUT 'Unspecified' annotations to allow peaks to continue pipeline.
Properly extracts FBgn IDs and looks up gene names from GTF.
"""

import argparse
import sys
import os
import subprocess
import tempfile
import re

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

def load_gene_names_from_gtf(gtf_file):
    """Load gene ID to gene name mapping from GTF file."""
    gene_mapping = {}

    if not gtf_file or not os.path.exists(gtf_file):
        print(f"Warning: GTF file {gtf_file} not found. Using gene IDs as gene names.")
        return gene_mapping

    print(f"Loading gene names from GTF: {gtf_file}")

    try:
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

                    feature = parts[2]
                    if feature == 'gene':
                        attributes = parse_gtf_attributes(parts[8])
                        gene_id = attributes.get('gene_id', '')
                        gene_name = attributes.get('gene_name', gene_id)

                        if gene_id:
                            gene_mapping[gene_id] = gene_name

                except Exception as e:
                    print(f"Warning: Error parsing GTF line {line_num}: {e}")
                    continue

        print(f"Loaded {len(gene_mapping)} gene name mappings from GTF")

    except Exception as e:
        print(f"Warning: Error reading GTF file: {e}")

    return gene_mapping

def extract_fbgn_from_crm_name(crm_name):
    """Extract FBgn ID from CRM name."""
    # Look for FBgn pattern in the CRM name
    fbgn_match = re.search(r'(FBgn\d+)', crm_name)
    if fbgn_match:
        return fbgn_match.group(1)

    # If no FBgn found, try other patterns
    if '_' in crm_name:
        # Try splitting and looking for gene-like identifiers
        parts = crm_name.split('_')
        for part in parts:
            if part.startswith('FBgn'):
                return part
            # Could add other gene ID patterns here if needed

    # If no recognizable gene ID, use the CRM name itself
    return crm_name

def run_bedtools_intersect(peaks_file, crm_file, output_file, overlap_fraction=0.0):
    """Run bedtools intersect and return results, filtering out Unspecified annotations."""
    try:
        # Create temporary file for bedtools output
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed') as temp_file:
            temp_output = temp_file.name

        # Run bedtools intersect with overlap fraction
        cmd = [
            'bedtools', 'intersect',
            '-a', peaks_file,
            '-b', crm_file,
            '-wa', '-wb'  # Write both A and B entries
        ]

        # Only add -f parameter if overlap_fraction > 0
        # Default bedtools behavior (no -f) = any overlap (≥1bp)
        if overlap_fraction > 0.0:
            cmd.extend(['-f', str(overlap_fraction)])

        print(f"Running: {' '.join(cmd)}")

        with open(temp_output, 'w') as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            print(f"Error running bedtools intersect: {result.stderr}")
            return False, []

        # Read and process results - FILTER OUT UNSPECIFIED
        intersections = []
        with open(temp_output, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 8:  # Peak (4 cols) + CRM (4+ cols)
                    # Peak information (first 4 columns)
                    peak_chr = parts[0]
                    peak_start = parts[1]
                    peak_end = parts[2]
                    peak_name = parts[3] if parts[3] != '.' else f"{peak_chr}:{peak_start}-{peak_end}"

                    # CRM information (remaining columns)
                    crm_chr = parts[4]
                    crm_start = parts[5]
                    crm_end = parts[6]
                    crm_name = parts[7] if len(parts) > 7 and parts[7] != '.' else f"{crm_chr}:{crm_start}-{crm_end}"

                    # FILTER OUT UNSPECIFIED ANNOTATIONS
                    if 'unspecified' in crm_name.lower():
                        print(f"Skipping unspecified CRM annotation: {crm_name}")
                        continue

                    intersections.append({
                        'peak_chr': peak_chr,
                        'peak_start': peak_start,
                        'peak_end': peak_end,
                        'peak_name': peak_name,
                        'crm_chr': crm_chr,
                        'crm_start': crm_start,
                        'crm_end': crm_end,
                        'crm_name': crm_name
                    })

        # Clean up temporary file
        os.unlink(temp_output)

        # Write intersections to output file
        with open(output_file, 'w') as f:
            f.write("peak_chr\tpeak_start\tpeak_end\tpeak_name\tcrm_chr\tcrm_start\tcrm_end\tcrm_name\n")
            for intersection in intersections:
                f.write(f"{intersection['peak_chr']}\t{intersection['peak_start']}\t{intersection['peak_end']}\t{intersection['peak_name']}\t")
                f.write(f"{intersection['crm_chr']}\t{intersection['crm_start']}\t{intersection['crm_end']}\t{intersection['crm_name']}\n")

        print(f"Found {len(intersections)} specified CRM intersections (filtered out unspecified)")
        return True, intersections

    except Exception as e:
        print(f"Error in bedtools intersect: {e}")
        return False, []

def annotate_crm_intersections(intersections, gene_mapping, output_file):
    """Annotate CRM intersections with gene information, properly extracting FBgn IDs."""
    try:
        # Create gene annotations from CRM regions
        gene_annotations = {}
        for intersection in intersections:
            crm_name = intersection['crm_name']
            peak_name = intersection['peak_name']

            # Extract FBgn ID from CRM name
            gene_id = extract_fbgn_from_crm_name(crm_name)

            # Get gene symbol from GTF mapping, fallback to gene_id
            gene_symbol = gene_mapping.get(gene_id, gene_id)

            if gene_id not in gene_annotations:
                gene_annotations[gene_id] = {
                    'gene_id': gene_id,
                    'gene_symbol': gene_symbol,
                    'annotation_type': 'crm',
                    'peak_names': [],
                    'crm_regions': []
                }

            gene_annotations[gene_id]['peak_names'].append(peak_name)
            gene_annotations[gene_id]['crm_regions'].append(crm_name)

        # Write annotated results
        with open(output_file, 'w') as f:
            f.write("gene_id\tgene_symbol\tannotation_type\tnum_peaks\tpeak_names\tcrm_regions\n")

            for gene_id, annotation in gene_annotations.items():
                gene_id_formatted = annotation['gene_id']
                gene_symbol = annotation['gene_symbol']
                annotation_type = annotation['annotation_type']
                num_peaks = len(set(annotation['peak_names']))
                peak_names = ','.join(sorted(set(annotation['peak_names'])))
                crm_regions = ','.join(sorted(set(annotation['crm_regions'])))

                f.write(f"{gene_id_formatted}\t{gene_symbol}\t{annotation_type}\t{num_peaks}\t{peak_names}\t{crm_regions}\n")

        print(f"Annotated {len(gene_annotations)} genes from specified CRM intersections")
        return True

    except Exception as e:
        print(f"Error in CRM annotation: {e}")
        return False

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Intersect peaks with CRM regions and annotate with gene information (excludes unspecified)'
    )
    parser.add_argument('--peaks', required=True, help='Input peaks file (BED format)')
    parser.add_argument('--crm', required=True, help='Input CRM regions file (BED format)')
    parser.add_argument('--gtf', help='GTF file for gene name lookup (optional)')
    parser.add_argument('--intersected', required=True, help='Output intersected peaks file')
    parser.add_argument('--annotated', required=True, help='Output annotated genes file')
    parser.add_argument('--overlap-fraction', type=float, default=0.0, help='Minimum overlap fraction (default: 0.0)')

    args = parser.parse_args()

    # Check input files
    if not os.path.exists(args.peaks):
        print(f"Error: Peaks file not found: {args.peaks}")
        sys.exit(1)

    if not os.path.exists(args.crm):
        print(f"Error: CRM file not found: {args.crm}")
        sys.exit(1)

    try:
        # Load gene names from GTF if provided
        gene_mapping = {}
        if args.gtf:
            gene_mapping = load_gene_names_from_gtf(args.gtf)

        print("Intersecting peaks with specified CRM regions...")
        success, intersections = run_bedtools_intersect(
            args.peaks, args.crm, args.intersected, args.overlap_fraction
        )

        if not success:
            print("Error: Failed to intersect peaks with CRM regions")
            sys.exit(1)

        if intersections:
            print("Annotating specified CRM intersections...")
            success = annotate_crm_intersections(intersections, gene_mapping, args.annotated)
            if not success:
                print("Error: Failed to annotate CRM intersections")
                sys.exit(1)
        else:
            # Create empty annotation file
            with open(args.annotated, 'w') as f:
                f.write("gene_id\tgene_symbol\tannotation_type\tnum_peaks\tpeak_names\tcrm_regions\n")
            print("No specified CRM intersections found")

        print("CRM intersection and annotation completed successfully!")
        print("Peaks intersecting only 'Unspecified' CRMs will continue to intron annotation.")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
