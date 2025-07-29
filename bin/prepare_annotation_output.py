#!/usr/bin/env python3
"""
Prepare final annotation output combining all annotation types with filtering and lncRNA-miRNA expansion.
- Creates cleaned raw output (raw minus problematic HOMER genes)
- Filters HOMER peaks that intersect with exons
- Removes targets associated with transposable elements, tRNA, rRNA, and pseudogenes
- Includes lncRNA-miRNA expansion
- Produces multiple output files with different filtering levels
"""

import argparse
import sys
import os
import re
from collections import defaultdict

def parse_gtf_file(gtf_file):
    """Parse GTF file to create gene_id to gene_symbol mapping and extract exon regions."""
    gene_mapping = {}
    gene_biotypes = {}
    exon_regions = set()

    if not os.path.exists(gtf_file):
        print(f"Warning: GTF file not found: {gtf_file}")
        return gene_mapping, gene_biotypes, exon_regions

    print(f"Parsing GTF file: {gtf_file}")

    try:
        with open(gtf_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue

                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 9:
                    continue


