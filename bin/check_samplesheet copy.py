#!/usr/bin/env python3

"""
Check samplesheet for multi-step peak annotation pipeline.
Updated to handle simplified samplesheet format without control column.
"""

import os
import sys
import errno
import argparse
import pandas as pd

def print_error(error, context='Line', context_str=''):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != '' and context_str != '':
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)

def check_samplesheet(file_in, file_out):
    """
    Check that the samplesheet follows the expected format.
    """

    # Check if file exists
    if not os.path.exists(file_in):
        print_error(f"Samplesheet file does not exist: {file_in}")

    # Read samplesheet
    try:
        df = pd.read_csv(file_in, dtype=str)
    except Exception as e:
        print_error(f"Failed to read samplesheet: {e}")

    # Check required columns
    required_columns = ['sample', 'replicate', 'peaks']

    for col in required_columns:
        if col not in df.columns:
            print_error(f"Missing required column: {col}")

    # Check for duplicate sample-replicate combinations
    if df.duplicated(subset=['sample', 'replicate']).any():
        print_error("Duplicate sample-replicate combinations found")

    # Check each row
    for index, row in df.iterrows():
        sample = row['sample']
        replicate = row['replicate']
        peaks = row['peaks']

        # Check sample name
        if pd.isna(sample) or sample == '':
            print_error("Sample name cannot be empty", f"Line {index + 2}")

        # Check replicate
        try:
            replicate_int = int(replicate)
            if replicate_int < 1:
                print_error("Replicate must be >= 1", f"Line {index + 2}")
        except (ValueError, TypeError):
            print_error("Replicate must be an integer", f"Line {index + 2}")

        # Check peaks file
        if pd.isna(peaks) or peaks == '':
            print_error("Peaks file path cannot be empty", f"Line {index + 2}")

        if not os.path.exists(peaks):
            print_error(f"Peaks file does not exist: {peaks}", f"Line {index + 2}")

        # Check if peaks file is in BED format (basic check)
        try:
            with open(peaks, 'r') as f:
                first_line = f.readline().strip()
                if first_line and not first_line.startswith('#'):
                    parts = first_line.split('\t')
                    if len(parts) < 3:
                        print_error(f"Peaks file does not appear to be in BED format (need at least 3 columns): {peaks}", f"Line {index + 2}")
                    # Check if coordinates are numeric
                    try:
                        int(parts[1])
                        int(parts[2])
                    except ValueError:
                        print_error(f"Invalid coordinates in peaks file: {peaks}", f"Line {index + 2}")
        except Exception as e:
            print_error(f"Failed to read peaks file {peaks}: {e}", f"Line {index + 2}")

    # Check that we have at least min_reps for each sample
    sample_counts = df['sample'].value_counts()
    for sample, count in sample_counts.items():
        if count < 2:
            print(f"WARNING: Sample '{sample}' has only {count} replicate(s). Consensus peaks require at least 2 replicates.")

    # Write validated samplesheet
    df.to_csv(file_out, index=False)
    print(f"Samplesheet validation passed. Output: {file_out}")

def main():
    parser = argparse.ArgumentParser(description='Check samplesheet for multi-step peak annotation pipeline')
    parser.add_argument('file_in', help='Input samplesheet')
    parser.add_argument('file_out', help='Output validated samplesheet')

    args = parser.parse_args()

    check_samplesheet(args.file_in, args.file_out)

if __name__ == "__main__":
    main()

# END OF SCRIPT
