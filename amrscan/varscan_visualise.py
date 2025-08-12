#!/usr/bin/env python3

# amrscan_pipeline/amrscan/varscan_visualise.py
import argparse
import os
import sys
import csv
import html
import re
from collections import defaultdict

class BColors:
    """A helper class to add color to terminal output."""
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    enabled = sys.stdout.isatty()

    @classmethod
    def _colorize(cls, color_code, text):
        return f"{color_code}{text}{cls.ENDC}" if cls.enabled else text
    @classmethod
    def cyan(cls, text): return cls._colorize(cls.OKCYAN, text)
    @classmethod
    def green(cls, text): return cls._colorize(cls.OKGREEN, text)
    @classmethod
    def red(cls, text): return cls._colorize(cls.FAIL, text)
    @classmethod
    def yellow(cls, text): return cls._colorize(cls.WARNING, text)

def sanitize_filename(name):
    """Converts a string into a safe filename."""
    name = name.replace(' ', '-').replace('/', '-')
    name = re.sub(r'[^\w\-.]', '', name)
    return name.lower() + ".html"

def load_metadata(metadata_path):
    """Loads metadata to map AROs to AMR_Gene_Families."""
    print(BColors.cyan(f"--- Loading metadata from: {metadata_path} ---"))
    if not os.path.exists(metadata_path):
        print(BColors.red(f"Error: Metadata file not found at '{metadata_path}'"), file=sys.stderr)
        sys.exit(1)
    
    aro_to_family = {}
    try:
        with open(metadata_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                aro_to_family[row['Sequence_ID']] = row['AMR_Gene_Family']
    except (Exception, KeyError) as e:
        print(BColors.red(f"Error reading or parsing metadata file: {e}"), file=sys.stderr)
        sys.exit(1)
    print(BColors.green(f"--- Loaded metadata for {len(aro_to_family)} AROs."))
    return aro_to_family

def load_variant_hits(hits_path):
    """Loads the _variant_hits.tsv file to get metadata for each read."""
    print(BColors.cyan(f"--- Loading variant hits from: {hits_path} ---"))
    if not os.path.exists(hits_path):
        print(BColors.red(f"Error: Variant hits file not found at '{hits_path}'"), file=sys.stderr)
        sys.exit(1)
    
    hits_by_query_id = {}
    try:
        with open(hits_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                hits_by_query_id[row['query_id']] = row
    except (Exception, KeyError) as e:
        print(BColors.red(f"Error reading or parsing variant hits file: {e}"), file=sys.stderr)
        sys.exit(1)
    
    print(BColors.green(f"--- Loaded metadata for {len(hits_by_query_id)} confirmed reads."))
    return hits_by_query_id

def load_variant_alignments(alignments_path):
    """Loads the pre-formatted alignment blocks from the _variant_alignments.txt file."""
    print(BColors.cyan(f"--- Loading alignments from: {alignments_path} ---"))
    if not os.path.exists(alignments_path):
        print(BColors.red(f"Error: Alignments file not found at '{alignments_path}'"), file=sys.stderr)
        sys.exit(1)

    alignments = {}
    current_query_id = None
    current_block_lines = []

    with open(alignments_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_query_id:
                    alignments[current_query_id] = "\n".join(current_block_lines)
                
                try:
                    current_query_id = line.split('|')[0].strip()[2:].strip()
                    current_block_lines = [line.strip()]
                except IndexError:
                    current_query_id = None
            elif current_query_id:
                current_block_lines.append(line.strip())
    
    if current_query_id:
        alignments[current_query_id] = "\n".join(current_block_lines)

    print(BColors.green(f"--- Loaded {len(alignments)} alignment blocks."))
    return alignments

def colorize_alignment_bases(ref_line, que_line):
    """
    Compares the REFSEQ and QUESEQ lines and returns an HTML string for the 
    query sequence with colored bases.
    """
    if not ref_line.startswith("REFSEQ:") or not que_line.startswith("QUESEQ:"):
        return html.escape(que_line)

    ref_seq = ref_line.split(":", 1)[1].strip()
    que_seq = que_line.split(":", 1)[1].strip()
    
    colored_que_html = ""
    for i in range(len(ref_seq)):
        ref_char = ref_seq[i]
        que_char = que_seq[i]

        if que_char == '.':
            colored_que_html += '.'
        elif ref_char == '-':
            colored_que_html += f'<span class="ins">{que_char}</span>'
        elif que_char == '-':
            colored_que_html += '<span class="del">-</span>'
        elif ref_char.upper() == que_char.upper():
            colored_que_html += f'<span class="match">{que_char}</span>'
        else:
            colored_que_html += f'<span class="mismatch">{que_char}</span>'
            
    return "QUESEQ:    " + colored_que_html

def generate_main_html_page(family_name, aro_html_blocks):
    """Generates the final, complete HTML page."""
    legend_html = """
    <div class="legend">
        <h3>Legend</h3>
        <div class="legend-item"><span class="match">A</span> Match</div>
        <div class="legend-item"><span class="mismatch">T</span> Mismatch</div>
        <div class="legend-item"><span class="ins">G</span> Insertion (in Query)</div>
        <div class="legend-item"><span class="del">-</span> Deletion (in Query)</div>
    </div>
    """

    return f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>AMRScan Variant Visualisation: {html.escape(family_name)}</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; margin: 2em; background-color: #f8f9fa; color: #212529; }}
        .container {{ max-width: 95%; margin: auto; }}
        h1 {{ color: #343a40; border-bottom: 2px solid #dee2e6; padding-bottom: 0.5em; }}
        h2 {{ color: #495057; margin-top: 2.5em; border-bottom: 1px solid #e9ecef; padding-bottom: 0.3em; }}
        .aro-section {{ margin-bottom: 3em; }}
        .alignment-block {{
            background-color: #fff;
            border-radius: 8px;
            padding: 1em;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
            margin-top: 1em;
            overflow-x: auto;
        }}
        pre {{
            font-family: "Menlo", "Consolas", "Courier New", monospace;
            font-size: 11px;
            line-height: 1.3;
            white-space: pre;
            margin: 0;
        }}
        .legend {{ background-color: #fff; padding: 1em; border-radius: 8px; margin: 2em 0; display: flex; gap: 2em; align-items: center; flex-wrap: wrap; }}
        .legend h3 {{ margin: 0; padding-right: 1em; border-right: 1px solid #dee2e6; }}
        .legend-item {{ display: flex; align-items: center; gap: 0.5em; font-size: 0.9em; }}
        .legend-item span {{ font-weight: bold; padding: 1px 3px; border-radius: 3px; }}
        .match    {{ color: #212529; }}
        .mismatch {{ color: #FFFFFF; background-color: #E74C3C; font-weight: bold; }}
        .ins      {{ color: #FFFFFF; background-color: #3498DB; font-weight: bold; }}
        .del      {{ color: #7f8c8d; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Variant Alignment Visualisation for {html.escape(family_name)}</h1>
        {legend_html}
        {''.join(aro_html_blocks)}
    </div>
</body>
</html>
"""

def main():
    parser = argparse.ArgumentParser(
        description="Generates HTML visualisations of confirmed variant alignments for all relevant gene families."
    )
    parser.add_argument(
        "--variant-hits",
        required=True,
        help="Path to the `_variant_hits.tsv` file from varscan_process_sam.py."
    )
    parser.add_argument(
        "--variant-alignments",
        required=True,
        help="Path to the `_variant_alignments.txt` file from varscan_process_sam.py."
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Path to the AMR database metadata file (TSV format)."
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Directory where all generated HTML visualisation files will be saved."
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Load all necessary data
    hits_data = load_variant_hits(args.variant_hits)
    alignment_data = load_variant_alignments(args.variant_alignments)
    aro_to_family = load_metadata(args.metadata)

    # 2. Group all hits by their gene family
    hits_by_family = defaultdict(list)
    for query_id, hit_info in hits_data.items():
        if query_id in alignment_data:
            aro_id = hit_info['ARO_matched']
            family_string = aro_to_family.get(aro_id, "Unknown_Family")
            
            # Combine hit metadata with its alignment block
            combined_info = hit_info.copy()
            combined_info['alignment_block'] = alignment_data[query_id]
            
            # Handle multi-family AROs
            for family in family_string.split(';'):
                hits_by_family[family].append(combined_info)
        else:
            print(BColors.yellow(f"Warning: Query ID '{query_id}' found in hits file but not in alignments file. Skipping."), file=sys.stderr)

    if not hits_by_family:
        print(BColors.red("Error: No matching hits found or no families could be determined. Cannot generate visualisations."), file=sys.stderr)
        sys.exit(1)
    
    print(BColors.cyan(f"\n--- Found {len(hits_by_family)} gene families with confirmed variants to visualise. ---"))
    created_files = []

    # 3. Main Loop: Generate an HTML file for each family
    for family_name, family_hits in sorted(hits_by_family.items()):
        print(BColors.cyan(f"\n--- Processing family: {family_name} ---"))
        
        # Group the family's hits by their specific ARO
        hits_by_aro = defaultdict(list)
        for hit in family_hits:
            hits_by_aro[hit['ARO_matched']].append(hit)

        all_aro_html_blocks = []
        for aro_id in sorted(hits_by_aro.keys()):
            aro_hits = hits_by_aro[aro_id]
            aro_hits.sort(key=lambda x: float(x.get('nucleotide_pid', 0)), reverse=True)

            html_content = [f'<section class="aro-section"><h2>{html.escape(aro_id)}</h2>']
            html_content.append(f'<p>Found {len(aro_hits)} confirmed read(s) for this variant.</p>')

            for hit in aro_hits:
                escaped_block = html.escape(hit['alignment_block'])
                lines = escaped_block.splitlines()
                
                header_line, ref_line, que_line, marker_line = lines[0], "", "", ""
                for line in lines[1:]:
                    if line.startswith("REFSEQ:"): ref_line = line
                    elif line.startswith("QUESEQ:"): que_line = line
                    elif line.startswith("MARKER:"): marker_line = line
                
                colored_que_html = colorize_alignment_bases(ref_line, que_line)

                html_content.append('<div class="alignment-block"><pre>')
                html_content.append(f'<span style="font-weight: bold;">{header_line}</span>\n')
                html_content.append(f"{ref_line}\n")
                html_content.append(f"{colored_que_html}\n")
                if marker_line:
                    html_content.append(f"{marker_line}\n")
                html_content.append('</pre></div>')

            html_content.append('</section>')
            all_aro_html_blocks.append("".join(html_content))

        # 4. Generate and write the final page for this family
        final_html = generate_main_html_page(family_name, all_aro_html_blocks)
        output_filename = sanitize_filename(family_name)
        output_path = os.path.join(args.output_dir, output_filename)
        
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(final_html)
            print(BColors.green(f"--- Successfully generated: {output_path} ---"))
            created_files.append(output_path)
        except IOError as e:
            print(BColors.red(f"Error writing to output file {output_path}: {e}"), file=sys.stderr)

    print(BColors.green("\n\n--- All Processing Complete ---"))
    print("The following output files were generated:")
    for f in sorted(created_files): print(f"- {f}")

if __name__ == "__main__":
    main()