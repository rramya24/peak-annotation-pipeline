#!/usr/bin/env python3

"""
Intersect peaks with intron regions and annotate with gene information.
Preserves peak names throughout the process.
"""

import argparse
import sys
import os
import subprocess
import tempfile

def run_bedtools_intersect(peaks_file, intron_file, output_file, overlap_fraction=0.0):
    """Run bedtools intersect and return results."""
    try:
        # Create temporary file for bedtools output
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed') as temp_file:
            temp_output = temp_file.name

        # Run bedtools intersect with overlap fraction
        cmd = [
            'bedtools', 'intersect',
            '-a', peaks_file,
            '-b', intron_file,
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
                if len(parts) >= 8:  # Peak (4 cols) + Intron (4+ cols)
                    # Peak information (first 4 columns)
                    peak_chr = parts[0]
                    peak_start = parts[1]
                    peak_end = parts[2]
                    peak_name = parts[3] if parts[3] != '.' else f"{peak_chr}:{peak_start}-{peak_end}"

                    # Intron information (remaining columns)
                    intron_chr = parts[4]
                    intron_start = parts[5]
                    intron_end = parts[6]
                    intron_name = parts[7] if len(parts) > 7 and parts[7] != '.' else f"{intron_chr}:{intron_start}-{intron_end}"

                    intersections.append({
                        'peak_chr': peak_chr,
                        'peak_start': peak_start,
                        'peak_end': peak_end,
                        'peak_name': peak_name,
                        'intron_chr': intron_chr,
                        'intron_start': intron_start,
                        'intron_end': intron_end,
                        'intron_name': intron_name
                    })

        # Clean up temporary file
        os.unlink(temp_output)

        # Write intersections to output file
        with open(output_file, 'w') as f:
            f.write("peak_chr\tpeak_start\tpeak_end\tpeak_name\tintron_chr\tintron_start\tintron_end\tintron_name\n")
            for intersection in intersections:
                f.write(f"{intersection['peak_chr']}\t{intersection['peak_start']}\t{intersection['peak_end']}\t{intersection['peak_name']}\t")
                f.write(f"{intersection['intron_chr']}\t{intersection['intron_start']}\t{intersection['intron_end']}\t{intersection['intron_name']}\n")

        print(f"Found {len(intersections)} intron intersections")
        return True

    except Exception as e:
        print(f"Error in bedtools intersect: {e}")
        return False

def annotate_intron_intersections(intersections_file, output_file):
    """Annotate intron intersections with gene information."""
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
                        'intron_chr': parts[4],
                        'intron_start': parts[5],
                        'intron_end': parts[6],
                        'intron_name': parts[7]
                    })

        # Create gene annotations from intron regions
        gene_annotations = {}
        for intersection in intersections:
            intron_name = intersection['intron_name']
            peak_name = intersection['peak_name']

            # Extract gene name from intron name
            # Intron names are usually formatted as gene_name_intron_1, gene_name_first_intron, etc.
            if '_intron_' in intron_name.lower():
                gene_name = intron_name.split('_intron_')[0]
            elif '_first_intron' in intron_name.lower():
                gene_name = intron_name.split('_first_intron')[0]
            elif '_' in intron_name:
                # Take the first part before underscore as gene name
                gene_name = intron_name.split('_')[0]
            else:
                gene_name = intron_name

            if gene_name not in gene_annotations:
                gene_annotations[gene_name] = {
                    'gene_id': f"INTRON_{gene_name}",
                    'gene_symbol': gene_name,
                    'annotation_type': 'intron',
                    'peak_names': [],
                    'intron_regions': []
                }

            gene_annotations[gene_name]['peak_names'].append(peak_name)
            gene_annotations[gene_name]['intron_regions'].append(intron_name)

        # Write annotated results
        with open(output_file, 'w') as f:
            f.write("gene_id\tgene_symbol\tannotation_type\tnum_peaks\tpeak_names\tintron_regions\n")

            for gene_name, annotation in gene_annotations.items():
                gene_id = annotation['gene_id']
                gene_symbol = annotation['gene_symbol']
                annotation_type = annotation['annotation_type']
                num_peaks = len(set(annotation['peak_names']))  # Remove duplicates
                peak_names = ','.join(sorted(set(annotation['peak_names'])))
                intron_regions = ','.join(sorted(set(annotation['intron_regions'])))

                f.write(f"{gene_id}\t{gene_symbol}\t{annotation_type}\t{num_peaks}\t{peak_names}\t{intron_regions}\n")

        print(f"Annotated {len(gene_annotations)} genes from intron intersections")
        print(f"Results written to: {output_file}")

        return True

    except Exception as e:
        print(f"Error in intron annotation: {e}")
        return False

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Intersect peaks with intron regions and annotate with gene information'
    )
    parser.add_argument(
        '--peaks',
        required=True,
        help='Input peaks file (BED format)'
    )
    parser.add_argument(
        '--introns',
        required=True,
        help='Input intron regions file (BED format)'
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

    if not os.path.exists(args.introns):
        print(f"Error: Intron file not found: {args.introns}")
        sys.exit(1)

    try:
        print("Intersecting peaks with intron regions...")
        success = run_bedtools_intersect(
            args.peaks,
            args.introns,
            args.intersected,
            args.overlap_fraction
        )

        if not success:
            print("Error: Failed to intersect peaks with intron regions")
            sys.exit(1)

        print("Annotating intron intersections...")
        success = annotate_intron_intersections(args.intersected, args.annotated)

        if not success:
            print("Error: Failed to annotate intron intersections")
            sys.exit(1)

        print("Intron intersection and annotation completed successfully!")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

# END OF SCRIPT
