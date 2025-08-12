#!/usr/bin/env python3

# amrscan_pipeline/amrscan/homscan_visualise.py
import argparse
import os
import sys
import csv
import html
import math
import json
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
    family_to_aros = defaultdict(list)
    try:
        with open(metadata_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                aro_id = row['Sequence_ID'] 
                family = row['AMR_Gene_Family'] 
                aro_to_family[aro_id] = family
                for fam in family.split(';'):
                    family_to_aros[fam].append(aro_id)
    except (Exception, KeyError) as e:
        print(BColors.red(f"Error reading or parsing metadata file: {e}"), file=sys.stderr)
        sys.exit(1)
    print(BColors.green(f"--- Loaded metadata for {len(family_to_aros)} families and {len(aro_to_family)} AROs."))
    return aro_to_family, family_to_aros

def load_reference_sequences(db_path, aro_ids):
    """Loads sequences for a specific list of AROs."""
    print(BColors.cyan(f"--- Loading reference sequences for {len(aro_ids)} AROs ---"))
    aro_sequences = {}
    aro_ids_set = set(aro_ids)
    
    with open(db_path, 'r') as f:
        current_id = None
        seq_parts = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id in aro_ids_set:
                    aro_sequences[current_id] = "".join(seq_parts)
                current_id = line[1:] 
                seq_parts = []
            elif current_id in aro_ids_set:
                seq_parts.append(line)
    
    if current_id in aro_ids_set:
        aro_sequences[current_id] = "".join(seq_parts)
        
    print(BColors.green(f"--- Found sequences for {len(aro_sequences)} of the requested AROs."))
    return aro_sequences

def load_all_hits(input_files, pid_cutoff, pid_type, aro_to_family):
    """
    Loads all hits from multiple files, filters by the chosen PID type, and enriches with family info.
    """
    pid_column = f"{pid_type}_pid"
    print(BColors.cyan(f"--- Loading all hits with {pid_column} >= {pid_cutoff*100:.2f}% ---"))
    all_hits_raw = []
    pid_cutoff_percent = pid_cutoff * 100.0
    for file_path in input_files:
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    try:
                        if float(row[pid_column]) >= pid_cutoff_percent:
                            row['AMR_Gene_Family'] = aro_to_family.get(row['ARO_matched'], 'Unknown_Family')
                            all_hits_raw.append(row)
                    except (ValueError, KeyError):
                        continue
        except FileNotFoundError:
            print(BColors.yellow(f"Warning: Input file not found, skipping: {file_path}"), file=sys.stderr)
    
    hits_by_query_id = defaultdict(list)
    for hit in all_hits_raw:
        simplified_hit = {
            'ARO_matched': hit['ARO_matched'],
            'AMR_Gene_Family': hit['AMR_Gene_Family'],
            'nucleotide_pid': hit['nucleotide_pid'],
            'protein_pid': hit['protein_pid'],
            'uniqueness': hit.get('uniqueness', 'UNKNOWN'),
            'position_on_ref': hit['position_on_ref'],
            'nucleotide_denominator': hit['nucleotide_denominator'],
        }
        hits_by_query_id[hit['query_id']].append(simplified_hit)

    print(BColors.green(f"--- Loaded {len(all_hits_raw)} hits from {len(hits_by_query_id)} unique reads passing the filter."))
    return all_hits_raw, hits_by_query_id

def load_wta_summary(wta_summary_path, target_gene_family):
    """Loads the WTA summary to identify AROs that were assigned reads and their counts for a specific family."""
    wta_winners = set()
    wta_counts = defaultdict(int)
    if not os.path.exists(wta_summary_path):
        return wta_winners, wta_counts

    try:
        with open(wta_summary_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                full_id_col = row['#AMR_Gene_Family;ARO']
                family_string_from_summary, aro_id = full_id_col.split(';', 1)
                if target_gene_family in set(family_string_from_summary.split(';')):
                    wta_winners.add(aro_id)
                    wta_counts[aro_id] = int(row['Count'])
    except (Exception, KeyError, ValueError) as e:
        print(BColors.red(f"Error reading or parsing WTA summary file: {e}"), file=sys.stderr)
        return set(), defaultdict(int)

    return wta_winners, wta_counts

def load_wta_assignments(wta_assignments_path):
    """Loads the WTA assignments file to map query_id to its winner ARO."""
    if not os.path.exists(wta_assignments_path):
        return {}
    
    query_id_to_wta_aro = {}
    try:
        with open(wta_assignments_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                query_id_to_wta_aro[row['query_id']] = row['winner_ARO']
    except (Exception, KeyError) as e:
        print(BColors.red(f"Error reading or parsing WTA assignments file: {e}"), file=sys.stderr)
        return {}
    return query_id_to_wta_aro

def get_color_for_pid(pid):
    """
    Returns a simplified, static color based on the PID value.
    Blue for 100%, Green for < 100%.
    """
    if pid >= 100.0:
        return "blue"
    else:
        return "#28a745"

def generate_single_aro_plot_html(aro_id, ref_seq_len, hits, is_wta_winner, original_read_count, wta_read_count, pid_type):
    """Generates the HTML for a single ARO's plot, including stats and SVG."""
    stats = {
        "ref_len": f"{ref_seq_len:,}",
        "num_reads_all": f"{original_read_count:,}",
        "num_reads_wta": f"{wta_read_count:,}",
        "coverage": "N/A"
    }

    hits_for_initial_layout = sorted(hits, key=lambda x: float(x['nucleotide_pid']), reverse=True)
    rect_height, gap = 8, 2
    lanes = [0]
    
    js_hits_data = []
    for hit in hits_for_initial_layout:
        start = int(hit['position_on_ref'])
        length = int(hit['nucleotide_denominator'])
        assigned_lane = -1
        for i, end_pos in enumerate(lanes):
            if start >= end_pos:
                assigned_lane = i
                break
        if assigned_lane == -1:
            assigned_lane = len(lanes)
            lanes.append(0)
        lanes[assigned_lane] = start + length
        hit_copy = hit.copy()
        hit_copy['draw_lane'] = assigned_lane
        js_hits_data.append(hit_copy)

    max_lanes = len(lanes) if lanes else 1
    svg_elements, axis_elements = [], []
    svg_filter_defs = """<defs><filter id="goldGlow" x="-50%" y="-50%" width="200%" height="200%"><feDropShadow dx="0" dy="0" stdDeviation="3" flood-color="#FFD700" flood-opacity="0.8"/></filter></defs>"""

    pid_column_for_color = f"{pid_type}_pid"

    for hit in js_hits_data:
        pid_for_color = float(hit[pid_column_for_color])
        color = get_color_for_pid(pid_for_color)
        
        y_pos = (max_lanes - hit['draw_lane'] - 1) * (rect_height + gap)
        uniqueness = hit.get('uniqueness', 'UNKNOWN')
        stroke_attr = 'stroke="#B8860B" stroke-width="1.5"' if uniqueness == 'UNIQUE' else ''
        filter_attr = 'filter="url(#goldGlow)"' if uniqueness == 'UNIQUE' else ''
        
        tooltip_text_short = (
            f"Read: {html.escape(hit['query_id'])} ({uniqueness})<br>"
            f"Nuc PID: {hit['nucleotide_pid']}%<br>"
            f"Prot PID: {hit['protein_pid']}%"
        )
        
        svg_elements.append(
            f'<rect x="{hit["position_on_ref"]}" y="{y_pos}" '
            f'width="{hit["nucleotide_denominator"]}" height="{rect_height}" '
            f'fill="{color}" rx="1" {stroke_attr} {filter_attr} '
            f'data-query-id="{html.escape(hit["query_id"])}"'
            f'data-aro-matched="{html.escape(hit["ARO_matched"])}"'
            f'data-wta-assigned-aro="{html.escape(hit.get("wta_assigned_aro", ""))}"'
            f'data-nuc-pid="{hit["nucleotide_pid"]}"'
            f'data-prot-pid="{hit["protein_pid"]}"'
            f'data-pos-on-ref="{hit["position_on_ref"]}"'
            f'data-aln-len="{hit["nucleotide_denominator"]}"'
            f'data-tooltip-short="{tooltip_text_short}" '
            f'onmousemove="showTooltip(event)" onmouseout="hideTooltip()" '
            f'onclick="showReadDetails(event, \'{html.escape(hit["query_id"])}\')"></rect>'
        )

    axis_y_pos = max_lanes * (rect_height + gap) + 10
    tick_height = 5
    tick_interval = max(1, int(10 ** math.floor(math.log10(ref_seq_len / 2))))
    if ref_seq_len / tick_interval > 10: tick_interval *= 2
    if ref_seq_len / tick_interval < 4: tick_interval /= 2
    tick_interval = max(1, int(tick_interval))

    axis_elements.append(f'<line class="axis-line" x1="0" y1="{axis_y_pos}" x2="{ref_seq_len}" y2="{axis_y_pos}"></line>')
    for i in range(0, ref_seq_len + 1, tick_interval):
        axis_elements.append(f'<line class="tick" x1="{i}" y1="{axis_y_pos}" x2="{i}" y2="{axis_y_pos + tick_height}"></line>')
        axis_elements.append(f'<text class="tick-label" x="{i}" y="{axis_y_pos + tick_height + 12}">{i}</text>')
    
    axis_label_y = axis_y_pos + tick_height + 28
    axis_elements.append(f'<text class="axis-label" x="{ref_seq_len / 2}" y="{axis_label_y}">Reference Position (bp)</text>')

    svg_width = ref_seq_len
    svg_height = axis_label_y + 10
    winner_attr = ' data-is-winner="true"' if is_wta_winner else ''

    return f"""
    <section class="aro-plot" id="plot-{html.escape(aro_id)}"{winner_attr} data-aro-id="{html.escape(aro_id)}">
        <h2>{html.escape(aro_id)}{' <span class="winner-badge">(WTA Winner)</span>' if is_wta_winner else ''}</h2>
        <div class="stats-grid">
            <div class="stat-card"><h3>Reference Length</h3><p>{stats['ref_len']} bp</p></div>
            <div class="stat-card"><h3>Reads Mapped</h3><p class="reads-mapped-count" data-original-count="{original_read_count}" data-wta-count="{wta_read_count}">{stats['num_reads_all']}</p></div>
            <div class="stat-card"><h3>Lateral Coverage</h3><p class="lateral-coverage-value">{stats['coverage']}</p></div>
        </div>
        <div class="chart-container">
            <svg width="{svg_width}" height="{svg_height}" data-ref-len="{ref_seq_len}">
                {svg_filter_defs}
                {''.join(svg_elements)}
                {''.join(axis_elements)}
            </svg>
        </div>
    </section>
    """

def generate_main_html_page(family_name, plot_htmls, hits_by_query_id_json, all_plot_data, wta_counts, has_wta_winners, has_wta_assignments, pid_type):
    """Generates the final, complete HTML page containing all plots, a static legend, and details box."""
    
    legend_html = f"""
    <div class="legend">
        <h3>PID: {pid_type.capitalize()}</h3>
        <div class="legend-item"><span class="color-box" style="background-color: blue;"></span>100.0%</div>
        <div class="legend-item"><span class="color-box" style="background-color: #28a745;"></span>< 100.0%</div>
        <div class="legend-item" style="margin-left: 2em;"><span class="color-box" style="background-color: transparent; border: 1.5px solid #B8860B; box-shadow: 0 0 5px #FFD700;"></span>Unique Read</div>
    </div>
    """

    script_js = f"""
    <script>
        const allReadHits = {hits_by_query_id_json};
        const allPlotData = {json.dumps(all_plot_data)};
        const wtaCounts = {json.dumps(wta_counts)};
        const hasWtaWinners = {json.dumps(has_wta_winners)};
        const hasWtaAssignments = {json.dumps(has_wta_assignments)};
        const tooltip = document.getElementById('tooltip');
        const readDetailsOverlay = document.getElementById('read-details-overlay');
        const detailsReadIdSpan = document.getElementById('details-read-id');
        const detailsContentDiv = document.getElementById('details-content');
        const body = document.body;
        const RECT_HEIGHT = 8, GAP = 2, TICK_HEIGHT = 5, TICK_LABEL_OFFSET = 12, AXIS_LABEL_OFFSET = 28;

        function showTooltip(evt) {{
            const target = evt.target;
            const tooltipText = target.getAttribute('data-tooltip-short');
            if (!tooltipText) return;

            tooltip.style.display = 'block';
            tooltip.innerHTML = tooltipText.replace(/\\n/g, '<br/>');

            const tooltipWidth = tooltip.offsetWidth;
            const tooltipHeight = tooltip.offsetHeight;
            const viewportWidth = window.innerWidth;
            const viewportHeight = window.innerHeight;

            let x = evt.clientX + 15;
            let y = evt.clientY + 15;

            if (x + tooltipWidth > viewportWidth) {{
                x = evt.clientX - tooltipWidth - 15;
            }}
            if (y + tooltipHeight > viewportHeight) {{
                y = evt.clientY - tooltipHeight - 15;
            }}
            
            tooltip.style.left = x + 'px';
            tooltip.style.top = y + 'px';
        }}

        function hideTooltip() {{
            tooltip.style.display = 'none';
        }}

        function showReadDetails(event, queryId) {{
            event.stopPropagation();
            const hitsForThisRead = allReadHits[queryId];
            if (!hitsForThisRead || hitsForThisRead.length === 0) {{
                detailsContentDiv.innerHTML = '<p>No detailed data available for this read.</p>';
                return;
            }}
            detailsReadIdSpan.textContent = queryId;
            let contentHtml = '<p>This read mapped to the following AROs:</p><table><thead><tr><th>Gene Family</th><th>ARO</th><th>Nuc PID</th><th>Protein PID</th><th>Uniqueness</th><th>Position</th><th>Length</th><th>WTA Winner</th></tr></thead><tbody>';
            hitsForThisRead.sort((a, b) => parseFloat(b.nucleotide_pid) - parseFloat(a.nucleotide_pid));
            hitsForThisRead.forEach(hit => {{
                const currentPlotH2Element = event.target.closest('.aro-plot').querySelector('h2');
                const currentPlotAroId = currentPlotH2Element.firstChild.nodeValue.trim();
                let rowClassAttr = (hit.ARO_matched === currentPlotAroId) ? 'highlight-row' : '';
                const wtaWinnerLabel = (hasWtaAssignments && hit.is_wta_winner_for_read) ? '<span class="wta-winner-label">YES</span>' : 'NO';
                contentHtml += `<tr class="${{rowClassAttr}}"><td>${{hit.AMR_Gene_Family}}</td><td>${{hit.ARO_matched}}</td><td>${{hit.nucleotide_pid}}%</td><td>${{hit.protein_pid}}%</td><td>${{hit.uniqueness}}</td><td>${{hit.position_on_ref}}</td><td>${{hit.nucleotide_denominator}}</td><td>${{wtaWinnerLabel}}</td></tr>`;
            }});
            contentHtml += '</tbody></table>';
            detailsContentDiv.innerHTML = contentHtml;
            readDetailsOverlay.style.visibility = 'visible';
            readDetailsOverlay.style.opacity = 1;
        }}

        function hideReadDetails() {{
            readDetailsOverlay.style.opacity = 0;
            readDetailsOverlay.addEventListener('transitionend', function handler() {{
                readDetailsOverlay.style.visibility = 'hidden';
                readDetailsOverlay.removeEventListener('transitionend', handler);
            }});
        }}

        readDetailsOverlay.addEventListener('click', function(event) {{
            if (event.target === readDetailsOverlay) {{
                hideReadDetails();
            }}
        }});

        function calculateLateralCoverage(hits, refLength) {{
            if (!hits || hits.length === 0 || refLength === 0) return '0.00%';
            const coverageSet = new Set();
            hits.forEach(hit => {{
                const start = parseInt(hit.getAttribute('data-pos-on-ref'));
                const length = parseInt(hit.getAttribute('data-aln-len'));
                for (let i = start; i < start + length; i++) {{
                    coverageSet.add(i);
                }}
            }});
            return ((coverageSet.size / refLength) * 100).toFixed(2) + '%';
        }}

        function applyWtaReadFilteringAndRelayout() {{
            const isWinnerOnlyMode = body.classList.contains('show-winner-only');
            document.querySelectorAll('.aro-plot').forEach(plot => {{
                const plotAroId = plot.getAttribute('data-aro-id');
                const isPlotWinner = plot.getAttribute('data-is-winner') === 'true';
                const refLength = parseInt(plot.querySelector('svg').getAttribute('data-ref-len'));
                let currentAroHits = allPlotData.filter(hit => hit.ARO_matched === plotAroId);
                let visibleHits = (isWinnerOnlyMode && hasWtaAssignments) ? (isPlotWinner ? currentAroHits.filter(hit => hit.wta_assigned_aro === plotAroId) : []) : currentAroHits;
                visibleHits.sort((a, b) => parseFloat(b.nucleotide_pid) - parseFloat(a.nucleotide_pid));
                const lanes = [0];
                let maxLanes = 0;
                visibleHits.forEach(hit => {{
                    const start = parseInt(hit.position_on_ref);
                    const length = parseInt(hit.nucleotide_denominator);
                    let assignedLane = -1;
                    for (let i = 0; i < lanes.length; i++) {{
                        if (start >= lanes[i]) {{
                            assignedLane = i;
                            break;
                        }}
                    }}
                    if (assignedLane === -1) {{
                        assignedLane = lanes.length;
                        lanes.push(0);
                    }}
                    lanes[assignedLane] = start + length;
                    hit.draw_lane = assignedLane;
                    maxLanes = Math.max(maxLanes, assignedLane + 1);
                }});
                const svg = plot.querySelector('svg');
                const rects = svg.querySelectorAll('rect');
                const rectMap = new Map();
                rects.forEach(rect => rectMap.set(rect.getAttribute('data-query-id'), rect));
                visibleHits.forEach(hit => {{
                    const rect = rectMap.get(hit.query_id);
                    if (rect) {{
                        const y_pos = (maxLanes - hit.draw_lane - 1) * (RECT_HEIGHT + GAP);
                        rect.setAttribute('y', y_pos);
                        rect.style.display = 'block';
                    }}
                }});
                rects.forEach(rect => {{
                    const queryId = rect.getAttribute('data-query-id');
                    const isVisible = visibleHits.some(hit => hit.query_id === queryId && hit.ARO_matched === rect.getAttribute('data-aro-matched'));
                    if (!isVisible) {{
                        rect.style.display = 'none';
                    }}
                }});
                const newSvgHeight = (maxLanes * (RECT_HEIGHT + GAP)) + 10 + TICK_HEIGHT + TICK_LABEL_OFFSET + AXIS_LABEL_OFFSET;
                svg.setAttribute('height', newSvgHeight);
                const newAxisYPos = (maxLanes * (RECT_HEIGHT + GAP)) + 10;
                svg.querySelector('.axis-line').setAttribute('y1', newAxisYPos);
                svg.querySelector('.axis-line').setAttribute('y2', newAxisYPos);
                svg.querySelectorAll('.tick').forEach(tick => {{
                    tick.setAttribute('y1', newAxisYPos);
                    tick.setAttribute('y2', newAxisYPos + TICK_HEIGHT);
                }});
                svg.querySelectorAll('.tick-label').forEach(label => {{
                    label.setAttribute('y', newAxisYPos + TICK_HEIGHT + TICK_LABEL_OFFSET);
                }});
                svg.querySelector('.axis-label').setAttribute('y', newAxisYPos + TICK_HEIGHT + AXIS_LABEL_OFFSET);
                const readsMappedP = plot.querySelector('.reads-mapped-count');
                readsMappedP.textContent = isWinnerOnlyMode ? (wtaCounts[plotAroId] ? wtaCounts[plotAroId].toLocaleString() : '0') : parseInt(readsMappedP.getAttribute('data-original-count')).toLocaleString();
                plot.querySelector('.lateral-coverage-value').textContent = calculateLateralCoverage(Array.from(rects).filter(r => r.style.display !== 'none'), refLength);
            }});
        }}
        
        const toggleWinnerBtn = document.getElementById('toggleWinnerBtn');
        if (toggleWinnerBtn && hasWtaWinners) {{
            toggleWinnerBtn.addEventListener('click', function() {{
                body.classList.toggle('show-winner-only');
                if (body.classList.contains('show-winner-only')) {{
                    toggleWinnerBtn.textContent = 'Show All AROs';
                    toggleWinnerBtn.classList.add('active');
                }} else {{
                    toggleWinnerBtn.textContent = 'Show WTA Winners Only';
                    toggleWinnerBtn.classList.remove('active');
                }}
                applyWtaReadFilteringAndRelayout();
            }});
        }}
        document.addEventListener('DOMContentLoaded', applyWtaReadFilteringAndRelayout);
    </script>
    """

    return f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>AMRScan Visualisation: {html.escape(family_name)}</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; margin: 2em; background-color: #f8f9fa; color: #212529; }}
        .container {{ max-width: 1200px; margin: auto; }}
        h1 {{ color: #343a40; border-bottom: 2px solid #dee2e6; padding-bottom: 0.5em; }}
        h2 {{ color: #495057; margin-top: 2em; border-bottom: 1px solid #e9ecef; padding-bottom: 0.3em; display: flex; align-items: center; gap: 0.5em;}}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1em; margin: 1em 0; }}
        .stat-card {{ background-color: #fff; border-radius: 8px; padding: 1em; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }}
        .stat-card h3 {{ margin: 0 0 0.5em 0; color: #6c757d; font-size: 1em; font-weight: 600; }}
        .stat-card p {{ margin: 0; color: #343a40; font-size: 1.75em; font-weight: 700; }}
        .chart-container {{ position: relative; width: 100%; margin-top: 1em; overflow-x: auto; padding: 1em; background-color: #fff; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }}
        #tooltip {{ position: fixed; display: none; padding: 8px; background-color: rgba(0, 0, 0, 0.85); color: white; border-radius: 4px; font-size: 12px; pointer-events: none; white-space: pre; z-index: 1000; }}
        .axis-line, .tick {{ stroke: #adb5bd; stroke-width: 2; }}
        .tick-label {{ font-size: 10px; fill: #6c757d; text-anchor: middle; }}
        .axis-label {{ font-size: 12px; font-weight: 600; fill: #495057; text-anchor: middle; }}
        .aro-plot {{ margin-top: 3em; }}
        .legend {{ background-color: #fff; padding: 1em; border-radius: 8px; margin: 2em 0; display: flex; gap: 2em; align-items: center; flex-wrap: wrap; }}
        .legend h3 {{ margin: 0; padding-right: 1em; border-right: 1px solid #dee2e6; }}
        .legend-item {{ display: flex; align-items: center; gap: 0.5em; font-size: 0.9em; }}
        .color-box {{ width: 15px; height: 15px; border-radius: 3px; border: 1px solid #ccc; }}
        .winner-badge {{ font-size: 0.7em; background-color: #28a745; color: white; padding: 0.2em 0.5em; border-radius: 4px; margin-left: 0.5em; }}
        #read-details-overlay {{ position: fixed; top: 0; left: 0; width: 100%; height: 100%; background-color: rgba(0, 0, 0, 0.5); display: flex; justify-content: center; align-items: center; z-index: 2000; visibility: hidden; opacity: 0; transition: opacity 0.3s ease; }}
        #read-details-box {{ background-color: white; padding: 25px; border-radius: 10px; box-shadow: 0 5px 15px rgba(0, 0, 0, 0.3); max-width: 1200px; max-height: 80vh; overflow-y: auto; position: relative; }}
        #read-details-box h3 {{ margin-top: 0; color: #343a40; border-bottom: 1px solid #dee2e6; padding-bottom: 10px; font-size: 1.5em; }}
        #read-details-box .close-btn {{ position: absolute; top: 10px; right: 15px; font-size: 1.5em; cursor: pointer; color: #6c757d; }}
        #read-details-box table {{ width: 100%; border-collapse: collapse; margin-top: 15px; }}
        #read-details-box th, #read-details-box td {{ border: 1px solid #e9ecef; padding: 8px 12px; text-align: left; }}
        #read-details-box th {{ background-color: #f1f3f5; font-weight: 600; }}
        #read-details-box tr:nth-child(even) {{ background-color: #f8f9fa; }}
        #read-details-box .highlight-row {{ background-color: #e6f7ff !important; font-weight: bold; }}
        #read-details-box .wta-winner-label {{ color: #28a745; font-weight: bold; margin-left: 5px; }}
        .controls button {{ background-color: #007bff; color: white; border: none; padding: 10px 20px; border-radius: 5px; cursor: pointer; font-size: 1em; transition: background-color 0.2s ease; }}
        .controls button:hover {{ background-color: #0056b3; }}
        .controls button.active {{ background-color: #28a745; }}
        .controls button:disabled {{ background-color: #6c757d; cursor: not-allowed; }}
        body.show-winner-only .aro-plot:not([data-is-winner="true"]) {{ display: none; }}
        .aro-plot rect.hidden-by-wta {{ display: none; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Coverage for Gene Family: {html.escape(family_name)}</h1>
        <!--
        <div class="controls">
            <button id="toggleWinnerBtn" {'disabled' if not has_wta_winners else ''}>{'Show WTA Winners Only' if has_wta_winners else 'No WTA Winners to Filter'}</button>
        </div>
        -->
        {legend_html}
        <div id="plots-container">{''.join(plot_htmls)}</div>
    </div>
    <div id="tooltip"></div>
    <div id="read-details-overlay">
        <div id="read-details-box">
            <span class="close-btn" onclick="hideReadDetails()">×</span>
            <h3>Read Details: <span id="details-read-id"></span></h3>
            <div id="details-content"></div>
        </div>
    </div>
    {script_js}
</body>
</html>
"""

def main():
    parser = argparse.ArgumentParser(
        description="Generates interactive HTML visualisations of read coverage for all relevant AMR gene families."
    )
    parser.add_argument("-i", "--input-files", required=True, help="Comma-delimited list of input `_hits.tsv` files.")
    parser.add_argument("-d", "--db", required=True, help="Path to the reference database FASTA file.")
    parser.add_argument("--metadata", required=True, help="Path to the AMR database metadata file (TSV format).")
    parser.add_argument("-o", "--output-dir", required=True, help="Directory where all generated HTML visualisation files will be saved.")
    parser.add_argument("--pid-cutoff", type=float, default=0.9, help="Minimum percent identity to consider a hit for visualisation (0.0-1.0 scale).")
    parser.add_argument("--pid-type", choices=['protein', 'nucleotide'], default='protein', help="PID type to use for filtering hits to display. Default: protein")
    parser.add_argument("--wta-summary", help="Optional: Path to the _summary_wta.tsv file. Enables filtering by WTA winners.")
    parser.add_argument("--wta-assignments", help="Optional: Path to the _wta_assignments.tsv file. Enables filtering of reads within plots.")
    args = parser.parse_args()

    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    # --- Data Loading ---
    aro_to_family, family_to_aros = load_metadata(args.metadata)
    all_hits_raw, all_hits_by_query_id = load_all_hits(input_files_list, args.pid_cutoff, args.pid_type, aro_to_family)
    
    if not all_hits_raw:
        print(BColors.yellow("Warning: No reads passed the PID cutoff. No visualisation files will be generated."))
        sys.exit(0)

    query_id_to_wta_aro = {}
    has_wta_assignments = False
    if args.wta_assignments:
        print(BColors.cyan(f"--- Loading WTA assignments from: {args.wta_assignments} ---"))
        query_id_to_wta_aro = load_wta_assignments(args.wta_assignments)
        if query_id_to_wta_aro:
            has_wta_assignments = True
            print(BColors.green(f"--- Loaded {len(query_id_to_wta_aro)} WTA assignments."))

    # --- Automatic Family Detection ---
    families_with_hits = set()
    for hit in all_hits_raw:
        for family in hit['AMR_Gene_Family'].split(';'):
            families_with_hits.add(family)
    
    print(BColors.cyan(f"\n--- Found {len(families_with_hits)} gene families with hits to visualise. ---"))

    created_files = []
    # --- Main Loop: Process and generate HTML for each family ---
    for family_name in sorted(list(families_with_hits)):
        print(BColors.cyan(f"\n--- Processing family: {family_name} ---"))

        target_aros = set(family_to_aros.get(family_name, []))
        family_hits_raw = [h for h in all_hits_raw if h['ARO_matched'] in target_aros]
        
        wta_winners, wta_counts = set(), defaultdict(int)
        has_wta_winners = False
        if args.wta_summary:
            wta_winners, wta_counts = load_wta_summary(args.wta_summary, family_name)
            if wta_winners: has_wta_winners = True
        
        family_hits_by_query_id = defaultdict(list)
        for hit in family_hits_raw:
            # *** FIX IS HERE ***
            # Provide a default empty string to .get() to prevent None values
            wta_winner_for_read = query_id_to_wta_aro.get(hit['query_id'], "")
            hit['wta_assigned_aro'] = wta_winner_for_read
            
            simplified_hit = all_hits_by_query_id[hit['query_id']]
            for s_hit in simplified_hit:
                s_hit['is_wta_winner_for_read'] = (s_hit['ARO_matched'] == wta_winner_for_read if wta_winner_for_read else False)
            family_hits_by_query_id[hit['query_id']] = simplified_hit

        hits_by_aro_for_plotting = defaultdict(list)
        for hit in family_hits_raw:
            hits_by_aro_for_plotting[hit['ARO_matched']].append(hit)

        aros_to_visualise = list(hits_by_aro_for_plotting.keys())
        if not aros_to_visualise: continue

        ref_sequences = load_reference_sequences(args.db, aros_to_visualise)
        
        plot_htmls = []
        all_plot_data_for_family = []
        for aro_id in sorted(aros_to_visualise):
            if aro_id not in ref_sequences: continue
            
            aro_hits = hits_by_aro_for_plotting[aro_id]
            ref_len = len(ref_sequences[aro_id])
            is_wta_winner = aro_id in wta_winners
            original_read_count = len(aro_hits)
            wta_read_count = wta_counts.get(aro_id, 0)

            plot_htmls.append(generate_single_aro_plot_html(aro_id, ref_len, aro_hits, is_wta_winner, original_read_count, wta_read_count, args.pid_type))
            
            for hit in aro_hits:
                all_plot_data_for_family.append({
                    'ARO_matched': hit['ARO_matched'], 'query_id': hit['query_id'],
                    'nucleotide_pid': float(hit['nucleotide_pid']), 'protein_pid': float(hit['protein_pid']),
                    'position_on_ref': int(hit['position_on_ref']),
                    'nucleotide_denominator': int(hit['nucleotide_denominator']), 'uniqueness': hit.get('uniqueness', 'UNKNOWN'),
                    'wta_assigned_aro': hit['wta_assigned_aro']
                })

        family_hits_by_query_id_json = json.dumps(family_hits_by_query_id)
        final_html = generate_main_html_page(family_name, plot_htmls, family_hits_by_query_id_json, all_plot_data_for_family, wta_counts, has_wta_winners, has_wta_assignments, args.pid_type)
        
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