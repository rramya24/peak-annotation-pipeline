#!/usr/bin/env python3

"""
Intersect peaks with CRM regions and annotate with gene information.
Preserves peak names throughout the process.
"""

import argparse
import sys
import os
import subprocess
import tempfile

def run_bedtools_intersect(peaks_file, crm_file, output_file, overlap_fraction=0.0):
    """Run bedtools intersect and return results."""
    try:
        # Create temporary file for bedtools output
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed') as temp_file:
            temp_output = temp_file.name

        # Run bedtools intersect with overlap fraction
        cmd = [
            'bedtools', 'intersect',
            '-a', peaks_file,
            '-b', crm_file,
            '-wa', '-wb',  # Write both A and B entries
            '-f', str(overlap_fraction)
        ]

        print(f"Running: {' '.join(cmd)}")

        with open(temp_output, 'w') as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            print(f"Error running bedtools intersect: {result.stderr}")
            return False

        # Read and process results
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

        print(f"Found {len(intersections)} CRM intersections")
        return True

    except Exception as e:
        print(f"Error in bedtools intersect: {e}")
        return False

def annotate_crm_intersections(intersections_file, output_file):
    """Annotate CRM intersections with gene information."""
    try:
        # Read intersections
        intersections = []
        with open(intersections_file, 'r') as f:
            header = f.readline().strip()
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 8:
                    intersections.append({
                        'peak_chr': parts[0],
                        'peak_start': parts[1],
                        'peak_end': parts[2],
                        'peak_name': parts[3],
                        'crm_chr': parts[4],
                        'crm_start': parts[5],
                        'crm_end': parts[6],
                        'crm_name': parts[7]
                    })

        # Create gene annotations from CRM regions
        # For CRM regions, we'll use the CRM name as the gene identifier
        gene_annotations = {}
        for intersection in intersections:
            crm_name = intersection['crm_name']
            peak_name = intersection['peak_name']

            # Extract gene name from CRM name if possible
            # CRM names might be formatted as gene_name_crm_1, gene_name_enhancer, etc.
            if '_crm_' in crm_name.lower():
                gene_name = crm_name.split('_crm_')[0]
            elif '_enhancer' in crm_name.lower():
                gene_name = crm_name.split('_enhancer')[0]
            elif '_' in crm_name:
                # Take the first part before underscore as gene name
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
                num_peaks = len(set(annotation['peak_names']))  # Remove duplicates
                peak_names = ','.join(sorted(set(annotation['peak_names'])))
                crm_regions = ','.join(sorted(set(annotation['crm_regions'])))

                f.write(f"{gene_id}\t{gene_symbol}\t{annotation_type}\t{num_peaks}\t{peak_names}\t{crm_regions}\n")

        print(f"Annotated {len(gene_annotations)} genes from CRM intersections")
        print(f"Results written to: {output_file}")

        return True

    except Exception as e:
        print(f"Error in CRM annotation: {e}")
        return False

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Intersect peaks with CRM regions and annotate with gene information'
    )
    parser.add_argument(
        '--peaks',
        required=True,
        help='Input peaks file (BED format)'
    )
    parser.add_argument(
        '--crm',
        required=True,
        help='Input CRM regions file (BED format)'
    )
    parser.add_argument(
        '--intersected',
        required=True,
        help='Output intersected peaks file'
    )
    parser.add_argument(
        '--annotated',
        required=True,
        help='Output annotated genes file'
    )
    parser.add_argument(
        '--overlap-fraction',
        type=float,
        default=0.0,
        help='Minimum overlap fraction (default: 0.0)'
    )

    args = parser.parse_args()

    # Check input files
    if not os.path.exists(args.peaks):
        print(f"Error: Peaks file not found: {args.peaks}")
        sys.exit(1)

    if not os.path.exists(args.crm):
        print(f"Error: CRM file not found: {args.crm}")
        sys.exit(1)

    try:
        print("Intersecting peaks with CRM regions...")
        success = run_bedtools_intersect(
            args.peaks,
            args.crm,
            args.intersected,
            args.overlap_fraction
        )

        if not success:
            print("Error: Failed to intersect peaks with CRM regions")
            sys.exit(1)

        print("Annotating CRM intersections...")
        success = annotate_crm_intersections(args.intersected, args.annotated)

        if not success:
            print("Error: Failed to annotate CRM intersections")
            sys.exit(1)

        print("CRM intersection and annotation completed successfully!")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

# END OF SCRIPT
