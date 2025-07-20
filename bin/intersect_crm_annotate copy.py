#!/usr/bin/env python3

"""
Intersect peaks with CRM regions and annotate with gene information.
FILTERS OUT 'Unspecified' annotations to allow peaks to continue pipeline.
"""

import argparse
import sys
import os
import subprocess
import tempfile

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

def get_peaks_with_specified_crm(peaks_file, intersections):
    """Get list of peak names that have specified CRM annotations."""
    specified_peaks = set()
    for intersection in intersections:
        specified_peaks.add(intersection['peak_name'])
    return specified_peaks

def create_non_intersected_peaks(peaks_file, specified_peaks, output_file):
    """Create file of peaks that did NOT intersect with specified CRMs."""
    try:
        with open(peaks_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 4:
                    peak_name = parts[3] if parts[3] != '.' else f"{parts[0]}:{parts[1]}-{parts[2]}"

                    # Only include peaks that did NOT intersect with specified CRMs
                    if peak_name not in specified_peaks:
                        outfile.write(line + '\n')

        print(f"Created non-intersected peaks file: {output_file}")
        return True

    except Exception as e:
        print(f"Error creating non-intersected peaks: {e}")
        return False

def annotate_crm_intersections(intersections, output_file):
    """Annotate CRM intersections with gene information."""
    try:
        # Create gene annotations from CRM regions
        gene_annotations = {}
        for intersection in intersections:
            crm_name = intersection['crm_name']
            peak_name = intersection['peak_name']

            # Extract gene name from CRM name if possible
            if '_' in crm_name and 'FBgn' in crm_name:
                # Handle FlyBase format: FBgn0034013_VT17121
                gene_name = crm_name.split('_')[0]
            elif '_crm_' in crm_name.lower():
                gene_name = crm_name.split('_crm_')[0]
            elif '_enhancer' in crm_name.lower():
                gene_name = crm_name.split('_enhancer')[0]
            elif '_' in crm_name:
                gene_name = crm_name.split('_')[0]
            else:
                gene_name = crm_name

            if gene_name not in gene_annotations:
                gene_annotations[gene_name] = {
                    'gene_id': f"CRM_{gene_name}",
                    'gene_symbol': gene_name,
                    'annotation_type': 'crm',
                    'peak_names': [],
                    'crm_regions': []
                }

            gene_annotations[gene_name]['peak_names'].append(peak_name)
            gene_annotations[gene_name]['crm_regions'].append(crm_name)

        # Write annotated results
        with open(output_file, 'w') as f:
            f.write("gene_id\tgene_symbol\tannotation_type\tnum_peaks\tpeak_names\tcrm_regions\n")

            for gene_name, annotation in gene_annotations.items():
                gene_id = annotation['gene_id']
                gene_symbol = annotation['gene_symbol']
                annotation_type = annotation['annotation_type']
                num_peaks = len(set(annotation['peak_names']))
                peak_names = ','.join(sorted(set(annotation['peak_names'])))
                crm_regions = ','.join(sorted(set(annotation['crm_regions'])))

                f.write(f"{gene_id}\t{gene_symbol}\t{annotation_type}\t{num_peaks}\t{peak_names}\t{crm_regions}\n")

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
        print("Intersecting peaks with specified CRM regions...")
        success, intersections = run_bedtools_intersect(
            args.peaks, args.crm, args.intersected, args.overlap_fraction
        )

        if not success:
            print("Error: Failed to intersect peaks with CRM regions")
            sys.exit(1)

        if intersections:
            print("Annotating specified CRM intersections...")
            success = annotate_crm_intersections(intersections, args.annotated)
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
