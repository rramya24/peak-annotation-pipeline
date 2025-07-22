#!/usr/bin/env python3

"""
Intersect peaks with intron regions and annotate with gene information.
Preserves peak names throughout the process.
Extracts FBgn IDs from full 4th column of the intron bed file and looks up gene names from GTF.
Preprocesses input files to use only first 4 columns.
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

def extract_fbgn_from_intron_name(intron_name):
    """
    Extract FBgn ID from intron name - for introns, the entire name IS the gene ID.
    Keep full FBgn format (no underscore splitting for introns).
    """
    # For intron files, the 4th column typically contains the full FBgn ID
    gene_id = intron_name.strip()
    return gene_id

def preprocess_bed_file(input_file, temp_file_prefix):
    """Create a temporary BED file with only the first 4 columns."""
    temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed', prefix=temp_file_prefix)

    try:
        with open(input_file, 'r') as infile, temp_file as outfile:
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                # Split by tabs and take only first 4 fields
                parts = line.split('\t')

                if len(parts) < 3:
                    print(f"Warning: Line {line_num} in {input_file} has only {len(parts)} fields, skipping")
                    continue

                # Take first 4 fields, pad with coordinate-based name if needed
                chr_col = parts[0]
                start_col = parts[1]
                end_col = parts[2]
                name_col = parts[3] if len(parts) >= 4 and parts[3].strip() != '' else f"{chr_col}:{start_col}-{end_col}"

                # Write 4-column BED format
                outfile.write(f"{chr_col}\t{start_col}\t{end_col}\t{name_col}\n")

        print(f"Preprocessed {input_file} -> {temp_file.name} (first 4 columns only)")
        return temp_file.name

    except Exception as e:
        print(f"Error preprocessing {input_file}: {e}")
        os.unlink(temp_file.name)
        return None

def run_bedtools_intersect_and_categorize(peaks_file, intron_file, intersected_output, non_intersected_output, overlap_fraction=0.0):
    """Run bedtools intersect and categorize peaks."""
    try:
        # Preprocess input files to ensure only first 4 columns are used
        print("Preprocessing input files to use first 4 columns only...")
        temp_peaks = preprocess_bed_file(peaks_file, "peaks_4col_")
        temp_introns = preprocess_bed_file(intron_file, "introns_4col_")

        if not temp_peaks or not temp_introns:
            print("Error: Failed to preprocess input files")
            return False, [], []

        # Create temporary file for bedtools output
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed') as temp_file:
            temp_output = temp_file.name

        # Run bedtools intersect with preprocessed files
        cmd = [
            'bedtools', 'intersect',
            '-a', temp_peaks,
            '-b', temp_introns,
            '-wa', '-wb'  # Write both A and B entries
        ]

        # Only add -f parameter if overlap_fraction > 0
        if overlap_fraction > 0.0:
            cmd.extend(['-f', str(overlap_fraction)])

        print(f"Running: bedtools intersect with preprocessed 4-column files")

        with open(temp_output, 'w') as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            print(f"Error running bedtools intersect: {result.stderr}")
            # Clean up temp files
            os.unlink(temp_peaks)
            os.unlink(temp_introns)
            return False, [], []

        # Read all peaks first
        all_peaks = {}
        with open(temp_peaks, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        peak_name = parts[3]
                        all_peaks[peak_name] = line

        print(f"Total input peaks: {len(all_peaks)}")

        # Read and process intersections
        intron_intersections = []
        peaks_with_intron = set()

        with open(temp_output, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) != 8:
                    continue

                # Extract fields - both files have 4 columns each
                peak_chr, peak_start, peak_end, peak_name = parts[0], parts[1], parts[2], parts[3]
                intron_chr, intron_start, intron_end, intron_name = parts[4], parts[5], parts[6], parts[7]

                # Show first few for verification
                if line_num <= 5:
                    print(f"Intersection {line_num}: Peak='{peak_name}', Intron='{intron_name}'")

                peaks_with_intron.add(peak_name)
                intron_intersections.append({
                    'peak_chr': peak_chr,
                    'peak_start': peak_start,
                    'peak_end': peak_end,
                    'peak_name': peak_name,
                    'intron_chr': intron_chr,
                    'intron_start': intron_start,
                    'intron_end': intron_end,
                    'intron_name': intron_name
                })

        # Create peaks for next step: non-intersected
        peaks_for_next_step = []

        for peak_name, peak_line in all_peaks.items():
            if peak_name not in peaks_with_intron:
                # This peak didn't intersect any introns
                peaks_for_next_step.append(peak_line)

        print(f"Peaks with intron intersections: {len(peaks_with_intron)}")
        print(f"Peaks for next step (non-intersected): {len(peaks_for_next_step)}")

        # Write intron intersections file (for annotation)
        with open(intersected_output, 'w') as f:
            f.write("peak_chr\tpeak_start\tpeak_end\tpeak_name\tintron_chr\tintron_start\tintron_end\tintron_name\n")
            for intersection in intron_intersections:
                f.write(f"{intersection['peak_chr']}\t{intersection['peak_start']}\t{intersection['peak_end']}\t{intersection['peak_name']}\t")
                f.write(f"{intersection['intron_chr']}\t{intersection['intron_start']}\t{intersection['intron_end']}\t{intersection['intron_name']}\n")

        # Write non-intersected peaks file (for next step)
        with open(non_intersected_output, 'w') as f:
            for peak_line in peaks_for_next_step:
                f.write(peak_line + '\n')

        # Clean up temporary files
        os.unlink(temp_peaks)
        os.unlink(temp_introns)
        os.unlink(temp_output)

        return True, intron_intersections, peaks_for_next_step

    except Exception as e:
        print(f"Error in bedtools intersect: {e}")
        return False, [], []

def annotate_intron_intersections(intersections, gene_mapping, output_file):
    """Annotate intron intersections with gene information."""
    try:
        # Create gene annotations from intron regions
        gene_annotations = {}
        for intersection in intersections:
            intron_name = intersection['intron_name']
            peak_name = intersection['peak_name']

            # Extract FBgn ID from intron name - keep full FBgn format
            gene_id = extract_fbgn_from_intron_name(intron_name)

            # Get gene symbol from GTF mapping, fallback to gene_id
            gene_symbol = gene_mapping.get(gene_id, gene_id)

            if gene_id not in gene_annotations:
                gene_annotations[gene_id] = {
                    'gene_id': gene_id,
                    'gene_symbol': gene_symbol,
                    'annotation_type': 'intron',
                    'peak_names': [],
                    'intron_regions': []
                }

            gene_annotations[gene_id]['peak_names'].append(peak_name)
            gene_annotations[gene_id]['intron_regions'].append(intron_name)

        # Write annotated results
        with open(output_file, 'w') as f:
            f.write("gene_id\tgene_symbol\tannotation_type\tnum_peaks\tpeak_names\tintron_regions\n")

            for gene_id, annotation in gene_annotations.items():
                gene_id_formatted = annotation['gene_id']
                gene_symbol = annotation['gene_symbol']
                annotation_type = annotation['annotation_type']
                num_peaks = len(set(annotation['peak_names']))  # Remove duplicates
                peak_names = ','.join(sorted(set(annotation['peak_names'])))
                intron_regions = ','.join(sorted(set(annotation['intron_regions'])))

                f.write(f"{gene_id_formatted}\t{gene_symbol}\t{annotation_type}\t{num_peaks}\t{peak_names}\t{intron_regions}\n")

        print(f"Annotated {len(gene_annotations)} genes from intron intersections")
        return True

    except Exception as e:
        print(f"Error in intron annotation: {e}")
        return False

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Intersect peaks with intron regions and annotate with gene information'
    )
    parser.add_argument('--peaks', required=True, help='Input peaks file (BED format)')
    parser.add_argument('--introns', required=True, help='Input intron regions file (BED format)')
    parser.add_argument('--gtf', help='GTF file for gene name lookup (optional)')
    parser.add_argument('--intersected', required=True, help='Output intersected peaks file')
    parser.add_argument('--annotated', required=True, help='Output annotated genes file')
    parser.add_argument('--non-intersected', required=True, help='Output non-intersected peaks file')
    parser.add_argument('--overlap-fraction', type=float, default=0.0, help='Minimum overlap fraction (default: 0.0)')

    args = parser.parse_args()

    # Check input files
    if not os.path.exists(args.peaks):
        print(f"Error: Peaks file not found: {args.peaks}")
        sys.exit(1)

    if not os.path.exists(args.introns):
        print(f"Error: Intron file not found: {args.introns}")
        sys.exit(1)

    try:
        # Load gene names from GTF if provided
        gene_mapping = {}
        if args.gtf:
            gene_mapping = load_gene_names_from_gtf(args.gtf)

        print("Intersecting peaks with intron regions...")
        success, intersections, peaks_for_next_step = run_bedtools_intersect_and_categorize(
            args.peaks, args.introns, args.intersected, args.non_intersected, args.overlap_fraction
        )

        if not success:
            print("Error: Failed to intersect peaks with intron regions")
            sys.exit(1)

        if intersections:
            print("Annotating intron intersections...")
            success = annotate_intron_intersections(intersections, gene_mapping, args.annotated)
            if not success:
                print("Error: Failed to annotate intron intersections")
                sys.exit(1)
        else:
            # Create empty annotation file
            with open(args.annotated, 'w') as f:
                f.write("gene_id\tgene_symbol\tannotation_type\tnum_peaks\tpeak_names\tintron_regions\n")
            print("No intron intersections found")

        print("Intron intersection and annotation completed successfully!")
        print(f"Peaks for next step: {len(peaks_for_next_step)} (non-intersected)")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
