#!/usr/bin/env python3

"""
Prepare final annotation output combining all annotation types.
Works with the existing PREPARE_ANNOTATION_OUTPUT module structure.
"""

import argparse
import sys
import os
from collections import defaultdict

def parse_homer_file(homer_file):
    """Parse HOMER annotation file and extract gene information."""
    genes = {}
    peaks_to_genes = {}

    if not os.path.exists(homer_file) or os.path.getsize(homer_file) == 0:
        print(f"Warning: HOMER file not found or empty: {homer_file}")
        return genes, peaks_to_genes

    try:
        with open(homer_file, 'r') as f:
            header = f.readline().strip()
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 16:  # HOMER format has many columns
                    peak_id = parts[0]
                    gene_name = parts[15] if len(parts) > 15 and parts[15] != '' else 'Unknown'
                    gene_type = parts[17] if len(parts) > 17 and parts[17] != '' else 'unknown'
                    annotation = parts[7] if len(parts) > 7 else 'Unknown'

                    if gene_name != 'Unknown' and gene_name != '':
                        if gene_name not in genes:
                            genes[gene_name] = {
                                'gene_id': gene_name,
                                'gene_symbol': gene_name,
                                'annotation_type': 'homer',
                                'peaks': [],
                                'gene_biotype': gene_type
                            }
                        genes[gene_name]['peaks'].append(peak_id)
                        peaks_to_genes[peak_id] = gene_name

    except Exception as e:
        print(f"Error parsing HOMER file: {e}")

    return genes, peaks_to_genes

def parse_bed_annotation_file(file_path, annotation_type):
    """Parse BED-style annotation files (CRM/intron)."""
    genes = {}
    peaks_to_genes = {}

    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        print(f"Warning: {annotation_type} file not found or empty: {file_path}")
        return genes, peaks_to_genes

    try:
        with open(file_path, 'r') as f:
            # Skip header if present
            header = f.readline().strip()
            if not header.startswith('#') and '\t' in header:
                # Process header line as data
                f.seek(0)

            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) >= 5:
                    gene_id = parts[0] if parts[0] != '' else 'Unknown'
                    gene_symbol = parts[1] if parts[1] != '' else gene_id
                    num_peaks = int(parts[3]) if parts[3].isdigit() else 0
                    peak_names = parts[4] if parts[4] != '' else ''
                    gene_biotype = parts[5] if len(parts) > 5 and parts[5] != '' else 'unknown'

                    if gene_symbol != 'Unknown' and gene_symbol != '':
                        if gene_symbol not in genes:
                            genes[gene_symbol] = {
                                'gene_id': gene_id,
                                'gene_symbol': gene_symbol,
                                'annotation_type': annotation_type,
                                'peaks': [],
                                'gene_biotype': gene_biotype
                            }

                        if peak_names:
                            peak_list = [p.strip() for p in peak_names.split(',') if p.strip()]
                            genes[gene_symbol]['peaks'].extend(peak_list)
                            for peak in peak_list:
                                peaks_to_genes[peak] = gene_symbol

    except Exception as e:
        print(f"Error parsing {annotation_type} file: {e}")

    return genes, peaks_to_genes

def combine_all_annotations(gene_symbol_files):
    """Combine annotations from all sources."""
    all_genes = {}
    all_peaks_to_genes = {}

    for file_path in gene_symbol_files:
        if not file_path or file_path == "NO_FILE_CRM" or file_path == "NO_FILE_INTRON" or file_path == "NO_FILE_HOMER":
            continue

        file_name = os.path.basename(file_path).lower()

        if 'homer' in file_name or 'annotatepeaks' in file_name:
            genes, peaks_to_genes = parse_homer_file(file_path)
            annotation_type = 'homer'
        elif 'crm' in file_name:
            genes, peaks_to_genes = parse_bed_annotation_file(file_path, 'crm')
            annotation_type = 'crm'
        elif 'intron' in file_name:
            genes, peaks_to_genes = parse_bed_annotation_file(file_path, 'intron')
            annotation_type = 'intron'
        else:
            # Try to parse as general annotation file
            genes, peaks_to_genes = parse_bed_annotation_file(file_path, 'unknown')
            annotation_type = 'unknown'

        print(f"Found {len(genes)} genes from {annotation_type} file: {file_path}")

        # Merge genes
        for gene_symbol, gene_data in genes.items():
            if gene_symbol not in all_genes:
                all_genes[gene_symbol] = {
                    'gene_id': gene_data['gene_id'],
                    'gene_symbol': gene_symbol,
                    'annotation_types': [annotation_type],
                    'peaks': gene_data['peaks'].copy(),
                    'gene_biotype': gene_data['gene_biotype']
                }
            else:
                # Merge with existing gene
                if annotation_type not in all_genes[gene_symbol]['annotation_types']:
                    all_genes[gene_symbol]['annotation_types'].append(annotation_type)
                all_genes[gene_symbol]['peaks'].extend(gene_data['peaks'])

        # Merge peak mappings
        all_peaks_to_genes.update(peaks_to_genes)

    # Remove duplicate peaks and sort
    for gene_symbol in all_genes:
        all_genes[gene_symbol]['peaks'] = sorted(list(set(all_genes[gene_symbol]['peaks'])))
        all_genes[gene_symbol]['num_peaks'] = len(all_genes[gene_symbol]['peaks'])
        all_genes[gene_symbol]['annotation_types'] = sorted(list(set(all_genes[gene_symbol]['annotation_types'])))

    return all_genes, all_peaks_to_genes

def write_output_files(all_genes, all_peaks_to_genes, output_prefix):
    """Write all output files."""

    # 1. All target genes file
    target_genes_file = f"{output_prefix}.all_target_genes.txt"
    with open(target_genes_file, 'w') as f:
        f.write("gene_id\tgene_symbol\tannotation_types\tnum_peaks\tpeak_names\tgene_biotype\n")
        for gene_symbol in sorted(all_genes.keys()):
            gene_data = all_genes[gene_symbol]
            annotation_types_str = ','.join(gene_data['annotation_types'])
            peak_names_str = ','.join(gene_data['peaks'])
            f.write(f"{gene_data['gene_id']}\t{gene_symbol}\t{annotation_types_str}\t")
            f.write(f"{gene_data['num_peaks']}\t{peak_names_str}\t{gene_data['gene_biotype']}\n")

    # 2. Peak annotation summary
    summary_file = f"{output_prefix}.peak_annotation_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Sample: {output_prefix}\n")
        f.write(f"Total target genes: {len(all_genes)}\n")
        f.write(f"Total annotated peaks: {len(all_peaks_to_genes)}\n")

        # Count by annotation type
        type_counts = defaultdict(int)
        for gene_data in all_genes.values():
            for ann_type in gene_data['annotation_types']:
                type_counts[ann_type] += 1

        f.write("\nGenes by annotation type:\n")
        for ann_type, count in sorted(type_counts.items()):
            f.write(f"  {ann_type}: {count}\n")

        # Count by biotype
        biotype_counts = defaultdict(int)
        for gene_data in all_genes.values():
            biotype_counts[gene_data['gene_biotype']] += 1

        f.write("\nGenes by biotype:\n")
        for biotype, count in sorted(biotype_counts.items()):
            f.write(f"  {biotype}: {count}\n")

    # 3. MultiQC file
    mqc_file = f"{output_prefix}.mqc.tsv"
    with open(mqc_file, 'w') as f:
        f.write("Sample\tTotal_Genes\tTotal_Peaks\tCRM_Genes\tIntron_Genes\tHOMER_Genes\n")

        type_counts = defaultdict(int)
        for gene_data in all_genes.values():
            for ann_type in gene_data['annotation_types']:
                type_counts[ann_type] += 1

        f.write(f"{output_prefix}\t{len(all_genes)}\t{len(all_peaks_to_genes)}\t")
        f.write(f"{type_counts.get('crm', 0)}\t{type_counts.get('intron', 0)}\t{type_counts.get('homer', 0)}\n")

    # 4. HTML report
    html_file = f"{output_prefix}.annotation_report.html"
    generate_html_report(all_genes, all_peaks_to_genes, html_file, output_prefix)

    return target_genes_file, summary_file, mqc_file, html_file

def generate_html_report(all_genes, all_peaks_to_genes, html_file, sample_name):
    """Generate HTML report."""
    total_genes = len(all_genes)
    total_peaks = len(all_peaks_to_genes)

    # Count by annotation type
    type_counts = defaultdict(int)
    biotype_counts = defaultdict(int)

    for gene_data in all_genes.values():
        for ann_type in gene_data['annotation_types']:
            type_counts[ann_type] += 1
        biotype_counts[gene_data['gene_biotype']] += 1

    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Peak Annotation Report - {sample_name}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
        .section {{ margin: 20px 0; }}
        .stats {{ display: flex; justify-content: space-around; margin: 20px 0; }}
        .stat-box {{ background-color: #e6f3ff; padding: 15px; border-radius: 5px; text-align: center; }}
        .stat-number {{ font-size: 24px; font-weight: bold; color: #0066cc; }}
        .stat-label {{ font-size: 14px; color: #666; }}
        table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
        .peak-name {{ font-family: monospace; background-color: #f8f8f8; padding: 2px 4px; border-radius: 3px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Peak Annotation Report</h1>
        <h2>Sample: {sample_name}</h2>
        <p>Multi-step peak annotation analysis results</p>
    </div>

    <div class="stats">
        <div class="stat-box">
            <div class="stat-number">{total_genes}</div>
            <div class="stat-label">Target Genes</div>
        </div>
        <div class="stat-box">
            <div class="stat-number">{total_peaks}</div>
            <div class="stat-label">Annotated Peaks</div>
        </div>
        <div class="stat-box">
            <div class="stat-number">{len(type_counts)}</div>
            <div class="stat-label">Annotation Types</div>
        </div>
    </div>

    <div class="section">
        <h3>Annotation Type Distribution</h3>
        <table>
            <tr><th>Annotation Type</th><th>Number of Genes</th></tr>
"""

    for ann_type, count in sorted(type_counts.items()):
        html_content += f"            <tr><td>{ann_type}</td><td>{count}</td></tr>\n"

    html_content += """
        </table>
    </div>

    <div class="section">
        <h3>Gene Biotype Distribution</h3>
        <table>
            <tr><th>Gene Biotype</th><th>Number of Genes</th></tr>
"""

    for biotype, count in sorted(biotype_counts.items()):
        html_content += f"            <tr><td>{biotype}</td><td>{count}</td></tr>\n"

    html_content += """
        </table>
    </div>

    <div class="section">
        <h3>Top 20 Target Genes (by number of peaks)</h3>
        <table>
            <tr><th>Gene Symbol</th><th>Number of Peaks</th><th>Peak Names (first 5)</th><th>Annotation Types</th><th>Gene Biotype</th></tr>
"""

    # Sort by number of peaks and show top 20
    sorted_genes = sorted(all_genes.items(), key=lambda x: x[1]['num_peaks'], reverse=True)[:20]

    for gene_symbol, gene_data in sorted_genes:
        peak_names_display = ', '.join(gene_data['peaks'][:5])
        if len(gene_data['peaks']) > 5:
            peak_names_display += ', ...'

        annotation_types_str = ', '.join(gene_data['annotation_types'])

        html_content += f"""            <tr>
                <td>{gene_symbol}</td>
                <td>{gene_data['num_peaks']}</td>
                <td><span class="peak-name">{peak_names_display}</span></td>
                <td>{annotation_types_str}</td>
                <td>{gene_data['gene_biotype']}</td>
            </tr>
"""

    html_content += """
        </table>
    </div>

</body>
</html>
"""

    with open(html_file, 'w') as f:
        f.write(html_content)

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Prepare final annotation output combining all annotation types'
    )
    parser.add_argument(
        '--gene_symbol_files',
        nargs='+',
        required=True,
        help='List of gene symbol annotation files (CRM, intron, HOMER)'
    )
    parser.add_argument(
        '--consensus_peaks',
        help='Consensus peaks file (not used in current implementation but kept for compatibility)'
    )
    parser.add_argument(
        '--output_prefix',
        required=True,
        help='Output file prefix'
    )

    args = parser.parse_args()

    try:
        print("Starting final annotation output preparation...")
        print(f"Input files: {args.gene_symbol_files}")
        print(f"Output prefix: {args.output_prefix}")

        # Combine all annotations
        all_genes, all_peaks_to_genes = combine_all_annotations(args.gene_symbol_files)

        if not all_genes:
            print("Warning: No genes found in any annotation files")
            # Create empty output files
            all_genes = {}
            all_peaks_to_genes = {}

        # Write output files
        target_genes_file, summary_file, mqc_file, html_file = write_output_files(
            all_genes, all_peaks_to_genes, args.output_prefix
        )

        print(f"\nFinal annotation output completed successfully!")
        print(f"Summary:")
        print(f"  Total target genes: {len(all_genes)}")
        print(f"  Total annotated peaks: {len(all_peaks_to_genes)}")
        print(f"  Output files:")
        print(f"    - Target genes: {target_genes_file}")
        print(f"    - Summary: {summary_file}")
        print(f"    - MultiQC: {mqc_file}")
        print(f"    - HTML report: {html_file}")

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
