#!/usr/bin/env python3

"""
Consensus peak calling script - Enhanced for multi-sample processing
Handles peak files with or without identifiers, generates proper names
"""

import os
import argparse
import pandas as pd
from collections import OrderedDict

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Find consensus peaks between replicates')
    parser.add_argument('--peak_files', nargs='+', help='Peak files for replicates')
    parser.add_argument('--min_reps', type=int, default=2, help='Minimum number of replicates')
    parser.add_argument('--prefix', default='consensus', help='Output prefix (sample name)')
    return parser.parse_args()

def read_peak_file_with_naming(peak_file, sample_name, replicate_num):
    """Read peak file and ensure proper naming"""
    peaks = []
    peak_counter = 1

    with open(peak_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])

                    # Handle peak naming
                    if len(parts) >= 4 and parts[3] != '.' and parts[3] != '':
                        # Use existing name
                        peak_name = parts[3]
                    else:
                        # Generate name based on sample and replicate
                        peak_name = f"{sample_name}_rep{replicate_num}_peak{peak_counter}"
                        peak_counter += 1

                    peaks.append((chrom, start, end, peak_name))

    return peaks

def main():
    args = parse_args()

    # Parse peak files with proper naming
    peak_files = args.peak_files
    min_reps = args.min_reps
    sample_name = args.prefix

    print(f"Processing {len(peak_files)} peak files for sample '{sample_name}'...")

    # Read peak files
    all_peaks = []
    for rep_idx, peak_file in enumerate(peak_files):
        replicate_num = rep_idx + 1
        peaks = read_peak_file_with_naming(peak_file, sample_name, replicate_num)

        for chrom, start, end, peak_name in peaks:
            all_peaks.append((chrom, start, end, peak_name, replicate_num))

        print(f"  Replicate {replicate_num}: {len(peaks)} peaks")

    # Sort by chromosome and position
    all_peaks.sort(key=lambda x: (x[0], x[1]))

    # Merge overlapping peaks and count replicates
    consensus_peaks = []
    current_interval = None
    current_replicates = set()
    peak_counter = 1

    for chrom, start, end, peak_name, replicate in all_peaks:
        if current_interval is None:
            current_interval = (chrom, start, end)
            current_replicates = {replicate}
        elif (chrom == current_interval[0] and
              start <= current_interval[2] + 150):  # 150bp merge distance
            # Merge intervals
            current_interval = (chrom, min(current_interval[1], start), max(current_interval[2], end))
            current_replicates.add(replicate)
        else:
            # Save previous interval if it meets criteria
            if len(current_replicates) >= min_reps:
                consensus_name = f"{sample_name}_consensus_peak_{peak_counter}"
                consensus_peaks.append((
                    current_interval[0],
                    current_interval[1],
                    current_interval[2],
                    consensus_name,
                    len(current_replicates)
                ))
                peak_counter += 1

            # Start new interval
            current_interval = (chrom, start, end)
            current_replicates = {replicate}

    # Don't forget the last interval
    if current_interval is not None and len(current_replicates) >= min_reps:
        consensus_name = f"{sample_name}_consensus_peak_{peak_counter}"
        consensus_peaks.append((
            current_interval[0],
            current_interval[1],
            current_interval[2],
            consensus_name,
            len(current_replicates)
        ))

    # Write consensus peaks
    output_file = f"{sample_name}.consensus_peaks.bed"
    with open(output_file, 'w') as f:
        f.write(f"# Consensus peaks for sample: {sample_name}\n")
        f.write(f"# chrom\tstart\tend\tname\tscore\tstrand\n")
        for chrom, start, end, name, score in consensus_peaks:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t.\n")

    # Write log
    log_file = f"{sample_name}.log"
    with open(log_file, 'w') as f:
        f.write(f"Consensus Peak Calling Summary for {sample_name}\n")
        f.write("="*50 + "\n")
        f.write(f"Sample name: {sample_name}\n")
        f.write(f"Input replicates: {len(peak_files)}\n")
        f.write(f"Minimum replicates required: {min_reps}\n")
        f.write(f"Consensus peaks found: {len(consensus_peaks)}\n")
        f.write("\n")

        # Count by replicate support
        rep_counts = {}
        for _, _, _, _, score in consensus_peaks:
            rep_counts[score] = rep_counts.get(score, 0) + 1

        f.write("Peaks by replicate support:\n")
        for rep_count in sorted(rep_counts.keys(), reverse=True):
            f.write(f"  {rep_count} replicates: {rep_counts[rep_count]} peaks\n")

        f.write("\n")
        f.write("Input files:\n")
        for i, peak_file in enumerate(peak_files):
            f.write(f"  Replicate {i+1}: {peak_file}\n")

    print(f"Consensus peaks for '{sample_name}': {len(consensus_peaks)} peaks")
    print(f"Output files: {output_file}, {log_file}")

if __name__ == "__main__":
    main()

# END OF SCRIPT
