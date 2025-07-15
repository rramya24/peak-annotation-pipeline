#!/usr/bin/env python3

"""
Prepare final annotation output combining all annotation types.
Preserves peak names and generates comprehensive reports.
"""

import argparse
import sys
import os
from collections import defaultdict

def parse_annotation_file(file_path, annotation_type):
    """Parse an annotation file and return gene annotations."""
    annotations = []

    if not os.path.exists(file_path):
        print(f"Warning: {annotation_type} file not found: {file_path}")
        return annotations

    try:
        with open(file_path, 'r') as f:
            header = f.readline().strip()
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 5:
                    gene_id = parts[0]
                    gene_symbol = parts[1]
                    file_annotation_type = parts[2]
                    num_peaks = int(parts[3]) if parts[3].isdigit() else 0
                    peak_names = parts[4]
                    gene_biotype = parts[5] if len(parts) > 5 else 'unknown'

                    annotations.append({
                        'gene_id': gene_id,
                        'gene_symbol': gene_symbol,
                        'annotation_type': annotation_type,
                        'num_peaks': num_peaks,
                        'peak_names': peak_names,
                        'gene_biotype': gene_biotype
                    })

    except Exception as e:
        print(f"Error parsing {annotation_type} file: {e}")

    return annotations

def combine_annotations(crm_file, intron_file, homer_file):
    """Combine all annotation types."""
    all_annotations = []

    # Parse each annotation type
    crm_annotations = parse_annotation_file(crm_file, 'crm')
    intron_annotations = parse_annotation_file(intron_file, 'intron')
    homer_annotations = parse_annotation_file(homer_file, 'homer')

    print(f"Found {len(crm_annotations)} CRM annotations")
    print(f"Found {len(intron_annotations)} intron annotations")
    print(f"Found {len(homer_annotations)} HOMER annotations")

    # Combine all annotations
    all_annotations.extend(crm_annotations)
    all_annotations.extend(intron_annotations)
    all_annotations.extend(homer_annotations)

    return all_annotations

def merge_gene_annotations(annotations):
    """Merge annotations for the same gene from different sources."""
    merged_genes = defaultdict(lambda: {
        'gene_id': '',
        'gene_symbol': '',
        'annotation_types': [],
        'total_peaks': 0,
        'peak_names': [],
        'gene_biotype': 'unknown'
    })

    for annotation in annotations:
        gene_symbol = annotation['gene_symbol']

        if not merged_genes[gene_symbol]['gene_id']:
            merged_genes[gene_symbol]['gene_id'] = annotation['gene_id']
            merged_genes[gene_symbol]['gene_symbol'] = gene_symbol
            merged_genes[gene_symbol]['gene_biotype'] = annotation['gene_biotype']

        # Add annotation type
        merged_genes[gene_symbol]['annotation_types'].append(annotation['annotation_type'])

        # Add peak information
        if annotation['peak_names']:
            peak_list = annotation['peak_names'].split(',')
            merged_genes[gene_symbol]['peak_names'].extend(peak_list)

        # Add to total peaks
        merged_genes[gene_symbol]['total_peaks'] += annotation['num_peaks']

    # Process merged data
    final_annotations = []
    for gene_symbol, data in merged_genes.items():
        # Remove duplicates and sort
        unique_annotation_types = sorted(set(data['annotation_types']))
        unique_peak_names = sorted(set(data['peak_names']))

        final_annotations.append({
            'gene_id': data['gene_id'],
            'gene_symbol': gene_symbol,
            'annotation_types': ','.join(unique_annotation_types),
            'num_peaks': len(unique_peak_names),
            'peak_names': ','.join(unique_peak_names),
            'gene_biotype': data['gene_biotype']
        })

    return final_annotations

def write_target_genes_file(annotations, output_file):
    """Write the target genes file."""
    with open(output_file, 'w') as f:
        f.write("gene_id\tgene_symbol\tannotation_types\tnum_peaks\tpeak_names\tgene_biotype\n")

        for annotation in annotations:
            f.write(f"{annotation['gene_id']}\t{annotation['gene_symbol']}\t{annotation['annotation_types']}\t")
            f.write(f"{annotation['num_peaks']}\t{annotation['peak_names']}\t{annotation['gene_biotype']}\n")

def write_peak_to_gene_mapping(annotations, output_file):
    """Write peak-to-gene mapping file."""
    peak_to_genes = defaultdict(list)

    # Build peak to gene mapping
    for annotation in annotations:
        if annotation['peak_names']:
            peak_list = annotation['peak_names'].split(',')
            for peak_name in peak_list:
                peak_name = peak_name.strip()
                if peak_name:
                    peak_to_genes[peak_name].append({
                        'gene_symbol': annotation['gene_symbol'],
                        'annotation_type': annotation['annotation_types'],
                        'gene_biotype': annotation['gene_biotype']
                    })

    # Write mapping
    with open(output_file, 'w') as f:
        f.write("peak_name\tgene_symbol\tannotation_types\tgene_biotype\n")

        for peak_name in sorted(peak_to_genes.keys()):
            genes = peak_to_genes[peak_name]
            for gene_info in genes:
                f.write(f"{peak_name}\t{gene_info['gene_symbol']}\t{gene_info['annotation_type']}\t{gene_info['gene_biotype']}\n")

def generate_html_report(annotations, peak_mapping_file, html_output, sample_name):
    """Generate HTML report with peak information."""
    # Calculate statistics
    total_genes = len(annotations)
    total_peaks = sum(len(set(a['peak_names'].split(','))) for a in annotations if a['peak_names'])

    # Count by annotation type
    annotation_counts = defaultdict(int)
    biotype_counts = defaultdict(int)

    for annotation in annotations:
        annotation_types = annotation['annotation_types'].split(',')
        for ann_type in annotation_types:
            annotation_counts[ann_type.strip()] += 1
        biotype_counts[annotation['gene_biotype']] += 1

    # Generate HTML
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
        <p>Generated on: {sample_name} peak annotation analysis</p>
    </div>

    <div class="stats">
        <div class="stat-box">
            <div class="stat-number">{total_genes}</div>
            <div class="stat-label">Target Genes</div>
        </div>
        <div class="stat-box">
            <div class="stat-number">{total_peaks}</div>
            <div class="stat-label">Total Peaks</div>
        </div>
        <div class="stat-box">
            <div class="stat-number">{len(annotation_counts)}</div>
            <div class="stat-label">Annotation Types</div>
        </div>
    </div>

    <div class="section">
        <h3>Annotation Type Distribution</h3>
        <table>
            <tr><th>Annotation Type</th><th>Number of Genes</th></tr>
"""

    for ann_type, count in sorted(annotation_counts.items()):
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
            <tr><th>Gene Symbol</th><th>Number of Peaks</th><th>Peak Names</th><th>Annotation Types</th><th>Gene Biotype</th></tr>
"""

    # Sort by number of peaks and show top 20
    sorted_annotations = sorted(annotations, key=lambda x: x['num_peaks'], reverse=True)[:20]

    for annotation in sorted_annotations:
        peak_names_html = []
        if annotation['peak_names']:
            peak_list = annotation['peak_names'].split(',')[:5]  # Show max 5 peaks
            for peak in peak_list:
                peak_names_html.append(f'<span class="peak-name">{peak.strip()}</span>')
            if len(annotation['peak_names'].split(',')) > 5:
                peak_names_html.append('...')

        html_content += f"""            <tr>
                <td>{annotation['gene_symbol']}</td>
                <td>{annotation['num_peaks']}</td>
                <td>{' '.join(peak_names_html)}</td>
                <td>{annotation['annotation_types']}</td>
                <td>{annotation['gene_biotype']}</td>
            </tr>
"""

    html_content += """
        </table>
    </div>

    <div class="section">
        <h3>All Target Genes</h3>
        <table>
            <tr><th>Gene Symbol</th><th>Gene ID</th><th>Peaks</th><th>Peak Names</th><th>Annotation Types</th><th>Biotype</th></tr>
"""

    for annotation in sorted(annotations, key=lambda x: x['gene_symbol']):
        peak_names_display = annotation['peak_names'][:100] + '...' if len(annotation['peak_names']) > 100 else annotation['peak_names']

        html_content += f"""            <tr>
                <td>{annotation['gene_symbol']}</td>
                <td>{annotation['gene_id']}</td>
                <td>{annotation['num_peaks']}</td>
                <td><span class="peak-name">{peak_names_display}</span></td>
                <td>{annotation['annotation_types']}</td>
                <td>{annotation['gene_biotype']}</td>
            </tr>
"""

    html_content += """
        </table>
    </div>

</body>
</html>
"""

    with open(html_output, 'w') as f:
        f.write(html_content)

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Prepare final annotation output combining all annotation types'
    )
    parser.add_argument(
        '--crm',
        help='CRM annotation file'
    )
    parser.add_argument(
        '--intron',
        help='Intron annotation file'
    )
    parser.add_argument(
        '--homer',
        help='HOMER annotation file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output target genes file'
    )
    parser.add_argument(
        '--peak-mapping',
        required=True,
        help='Output peak-to-gene mapping file'
    )
    parser.add_argument(
        '--html-report',
        required=True,
        help='Output HTML report'
    )
    parser.add_argument(
        '--sample-name',
        required=True,
        help='Sample name for reporting'
    )

    args = parser.parse_args()

    # Check that at least one annotation file is provided
    if not any([args.crm, args.intron, args.homer]):
        print("Error: At least one annotation file (--crm, --intron, or --homer) must be provided")
        sys.exit(1)

    try:
        print("Combining annotation files...")
        annotations = combine_annotations(
            args.crm or '',
            args.intron or '',
            args.homer or ''
        )

        if not annotations:
            print("Error: No annotations found in input files")
            sys.exit(1)

        print("Merging gene annotations...")
        merged_annotations = merge_gene_annotations(annotations)

        print("Writing target genes file...")
        write_target_genes_file(merged_annotations, args.output)

        print("Writing peak-to-gene mapping...")
        write_peak_to_gene_mapping(merged_annotations, args.peak_mapping)

        print("Generating HTML report...")
        generate_html_report(merged_annotations, args.peak_mapping, args.html_report, args.sample_name)

        print("Final annotation output completed successfully!")

        # Print summary
        print(f"\nSummary:")
        print(f"  Total target genes: {len(merged_annotations)}")
        print(f"  Total peaks: {sum(a['num_peaks'] for a in merged_annotations)}")
        print(f"  Output files:")
        print(f"    - Target genes: {args.output}")
        print(f"    - Peak mapping: {args.peak_mapping}")
        print(f"    - HTML report: {args.html_report}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

# END OF SCRIPT
