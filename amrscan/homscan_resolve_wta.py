#!/usr/bin/env python3

import argparse
import os
import sys
import csv
from collections import defaultdict, Counter

class BColors:
    """A helper class to add color to terminal output."""
    HEADER = '\033[95m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
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


def load_metadata(metadata_path):
    """Loads the metadata file to map ARO to AMR_Gene_Family."""
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
    print(BColors.green(f"--- Loaded {len(aro_to_family)} metadata entries."))
    return aro_to_family

def load_and_filter_hits(input_files, pid_cutoff, pid_type):
    """Loads all hits from multiple files and applies the chosen PID cutoff."""
    pid_column = f"{pid_type}_pid"
    print(BColors.cyan(f"--- Loading hits from {len(input_files)} file(s) with {pid_column} >= {pid_cutoff*100:.2f}% ---"))
    all_hits = []
    pid_cutoff_percent = pid_cutoff * 100.0
    for file_path in input_files:
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    try:
                        if float(row[pid_column]) >= pid_cutoff_percent:
                            all_hits.append(row)
                    except (ValueError, KeyError):
                        continue
        except FileNotFoundError:
            print(BColors.yellow(f"Warning: Input file not found, skipping: {file_path}"), file=sys.stderr)
    print(BColors.green(f"--- Loaded a total of {len(all_hits)} hits passing the filter."))
    return all_hits

def calculate_wta_metrics(hits):
    """Calculates lateral coverage and read counts for all AROs to create a global ranking."""
    print(BColors.cyan("--- Calculating lateral coverage and read counts for all AROs ---"))
    coverage_sets = defaultdict(set)
    read_counts = Counter()
    target_lengths = {}

    for hit in hits:
        aro = hit['ARO_matched']
        read_counts[aro] += 1
        try:
            pos = int(hit['position_on_ref'])
            aln_len = int(hit['nucleotide_denominator'])
            target_len = int(hit['target_length'])
            
            if aro not in target_lengths:
                target_lengths[aro] = target_len
            
            coverage_sets[aro].update(range(pos, pos + aln_len))
        except (ValueError, KeyError):
            continue

    aro_metrics = {}
    for aro, positions in coverage_sets.items():
        total_len = target_lengths.get(aro)
        if total_len and total_len > 0:
            coverage = (len(positions) / total_len) * 100.0
            aro_metrics[aro] = {'coverage': coverage, 'reads': read_counts[aro]}
    
    print(BColors.green(f"--- Calculated metrics for {len(aro_metrics)} AROs."))
    return aro_metrics

def get_local_winner(query_hits, all_aro_metrics, pid_type, debug_mode=False):
    """
    Determines the winner from a specific set of candidate AROs for a single query read.
    Tie-breaking order:
    1. Lateral Coverage (global metric, higher is better)
    2. Specific PID (protein or nucleotide, higher is better)
    3. Read Count (global metric, higher is better)
    4. Alphabetical ARO ID (for deterministic tie-breaking)
    """
    if not query_hits:
        return None

    pid_column = f"{pid_type}_pid"
    hits_by_candidate_aro = {hit['ARO_matched']: hit for hit in query_hits}
    candidate_aros = list(hits_by_candidate_aro.keys())
    sortable_candidates = []
    
    if debug_mode and query_hits:
        print(BColors.yellow(f"\n--- Debugging WTA resolution for query_id: {query_hits[0]['query_id']} ---"))

    for aro in candidate_aros:
        global_metrics = all_aro_metrics.get(aro)
        global_coverage = global_metrics['coverage'] if global_metrics else -1.0
        global_reads = global_metrics['reads'] if global_metrics else -1
        
        try:
            specific_pid = float(hits_by_candidate_aro[aro][pid_column])
        except (ValueError, KeyError):
            specific_pid = -1.0

        sort_key = (-global_coverage, -specific_pid, -global_reads, aro)
        sortable_candidates.append((sort_key, aro))

        if debug_mode:
            print(BColors.yellow(f"  Candidate ARO: {aro}"))
            print(BColors.yellow(f"    Global Coverage: {global_coverage:.8f}%"))
            print(BColors.yellow(f"    Specific {pid_type.capitalize()} PID: {specific_pid:.8f}%"))
            print(BColors.yellow(f"    Global Reads: {global_reads}"))
            print(BColors.yellow(f"    Sort Key: {sort_key}"))

    sortable_candidates.sort()
    winner_aro = sortable_candidates[0][1] if sortable_candidates else None
    
    if debug_mode:
        print(BColors.yellow(f"  Final Winner for this read: {winner_aro}"))
        if query_hits:
            print(BColors.yellow(f"--- End Debugging for query_id: {query_hits[0]['query_id']} ---\n"))

    return winner_aro

def resolve_and_count_hits_wta(hits, all_aro_metrics, pid_type, debug_mode=False):
    """
    Resolves every ambiguous hit using a local WTA rule based on global metrics.
    Returns final counts and a mapping of query_id to its assigned winner ARO.
    """
    print(BColors.cyan("--- Resolving all reads with 'Winner Takes All' logic ---"))
    hits_by_query = defaultdict(list)
    for hit in hits:
        hits_by_query[hit['query_id']].append(hit)
    
    print(BColors.green(f"--- Processing {len(hits_by_query)} unique reads in total. ---"))

    final_counts = Counter()
    wta_assignments = {}

    for query_id, query_hits in hits_by_query.items():
        candidate_aros = {h['ARO_matched'] for h in query_hits}
        
        if len(candidate_aros) == 1:
            winner = list(candidate_aros)[0]
        else:
            winner = get_local_winner(query_hits, all_aro_metrics, pid_type, debug_mode)
        
        if winner:
            final_counts[winner] += 1
            wta_assignments[query_id] = winner
    
    print(BColors.green("--- Read resolution complete. No 'multiple' categories remaining."))
    return final_counts, wta_assignments

def write_wta_summary(final_counts, aro_to_family, output_path):
    """Writes the final, resolved summary file, sorted alphabetically."""
    print(BColors.cyan(f"--- Writing final summary to: {output_path} ---"))
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["#AMR_Gene_Family;ARO", "Count"])
        
        summary_data = []
        for aro, count in final_counts.items():
            family = aro_to_family.get(aro, "Unknown_Family")
            summary_data.append((f"{family};{aro}", count))
            
        summary_data.sort(key=lambda x: x[0])
        writer.writerows(summary_data)
    print(BColors.green("--- Write complete."))

def write_wta_assignments(wta_assignments, output_path):
    """Writes the query_id to winner ARO assignments."""
    print(BColors.cyan(f"--- Writing WTA assignments to: {output_path} ---"))
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["query_id", "winner_ARO"])
        for query_id, winner_aro in wta_assignments.items():
            writer.writerow([query_id, winner_aro])
    print(BColors.green("--- WTA assignments write complete."))

def main():
    parser = argparse.ArgumentParser(
        description="Resolves all ambiguous AMR hits using a local 'Winner Takes All' (WTA) method."
    )
    parser.add_argument("-i", "--input-files", required=True, help="Comma-delimited list of input `_hits.tsv` files.")
    parser.add_argument("--metadata", required=True, help="Path to the AMR database metadata file (TSV format).")
    parser.add_argument("--tmp-dir", default=".", help="Directory to store the output report file.")
    parser.add_argument("--output-prefix", required=True, help="Prefix for the combined output summary file.")
    parser.add_argument("--pid-cutoff", type=float, default=0.9, help="Minimum percent identity to consider a hit (0.0-1.0 scale).")
    # <<< NEW ARGUMENT >>>
    parser.add_argument("--pid-type", choices=['protein', 'nucleotide'], default='protein', help="PID type to use for filtering and tie-breaking. Default: protein")
    parser.add_argument("--debug-wta", action="store_true", help="Enable verbose debugging output for WTA resolution logic.")
    args = parser.parse_args()

    if not 0.0 <= args.pid_cutoff <= 1.0:
        print(BColors.red("Error: --pid-cutoff must be a fraction between 0.0 and 1.0."), file=sys.stderr)
        sys.exit(1)
        
    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.tmp_dir, exist_ok=True)

    output_path_summary = os.path.join(args.tmp_dir, f"{args.output_prefix}_summary_wta.tsv")
    output_path_assignments = os.path.join(args.tmp_dir, f"{args.output_prefix}_wta_assignments.tsv")

    aro_to_family = load_metadata(args.metadata)
    all_hits = load_and_filter_hits(input_files_list, args.pid_cutoff, args.pid_type)

    if not all_hits:
        print(BColors.yellow("Warning: No hits passed the PID cutoff. Cannot generate a summary."))
        with open(output_path_summary, 'w') as f: f.write("#AMR_Gene_Family;ARO\tCount\n")
        with open(output_path_assignments, 'w') as f: f.write("query_id\twinner_ARO\n")
        sys.exit(0)

    all_aro_metrics = calculate_wta_metrics(all_hits)
    final_counts, wta_assignments = resolve_and_count_hits_wta(all_hits, all_aro_metrics, args.pid_type, args.debug_wta)

    write_wta_summary(final_counts, aro_to_family, output_path_summary)
    write_wta_assignments(wta_assignments, output_path_assignments)

    print(BColors.green("\n\n--- WTA Resolution Complete ---"))
    print(f"The following output files were generated:")
    print(f"- {output_path_summary}")
    print(f"- {output_path_assignments}")

if __name__ == "__main__":
    main()