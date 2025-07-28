#!/usr/bin/env python3

import argparse
import os
import sys
import csv
import re
from collections import defaultdict, Counter

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

def load_metadata(metadata_path):
    """Loads metadata to map Sequence_ID to AMR_Gene_Family and SeqNucLength."""
    print(BColors.cyan(f"--- Loading metadata from: {metadata_path} ---"))
    if not os.path.exists(metadata_path):
        print(BColors.red(f"Error: Metadata file not found at '{metadata_path}'"), file=sys.stderr)
        sys.exit(1)
    aro_to_family, aro_to_length = {}, {}
    try:
        with open(metadata_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                aro_to_family[row['Sequence_ID']] = row['AMR_Gene_Family']
                aro_to_length[row['Sequence_ID']] = int(row['SeqNucLength']) if row.get('SeqNucLength') else 0
    except (Exception, KeyError) as e:
        print(BColors.red(f"Error reading or parsing metadata file: {e}"), file=sys.stderr)
        sys.exit(1)
    print(BColors.green(f"--- Loaded metadata for {len(aro_to_family)} AROs."))
    return aro_to_family, aro_to_length

def load_uscg_metrics(uscg_report_path):
    """Loads the Overall RPK and FPK values from the USCG report."""
    print(BColors.cyan(f"--- Loading USCG metrics from: {uscg_report_path} ---"))
    if not os.path.exists(uscg_report_path):
        print(BColors.red(f"Error: USCG report file not found at '{uscg_report_path}'"), file=sys.stderr)
        return None, None
        
    uscg_rpk, uscg_fpk = None, None
    try:
        with open(uscg_report_path, 'r') as f:
            for line in f:
                if line.startswith("Overall_RPK_Across_All_USCGs"):
                    uscg_rpk = float(line.strip().split('\t')[1])
                elif line.startswith("Overall_FPK_Across_All_USCGs"):
                    uscg_fpk = float(line.strip().split('\t')[1])
    except (IOError, ValueError, IndexError) as e:
        print(BColors.red(f"Error reading or parsing USCG report file: {e}"), file=sys.stderr)
        return None, None
    
    if uscg_rpk is None:
        print(BColors.yellow("Warning: 'Overall_RPK_Across_All_USCGs' not found in the USCG report."), file=sys.stderr)
    if uscg_fpk is None:
        print(BColors.yellow("Warning: 'Overall_FPK_Across_All_USCGs' not found in the USCG report. FPK-based normalization will be skipped."), file=sys.stderr)
        
    return uscg_rpk, uscg_fpk

def load_total_bases(filepath):
    """Loads the total base count from a file."""
    print(BColors.cyan(f"--- Loading total bases from: {filepath} ---"))
    if not os.path.exists(filepath):
        print(BColors.red(f"Error: Total bases file not found at '{filepath}'"), file=sys.stderr)
        sys.exit(1)
    try:
        with open(filepath, 'r') as f:
            bases = int(f.read().strip())
            print(BColors.green(f"--- Found total bases: {bases:,} ---"))
            return bases
    except (IOError, ValueError) as e:
        print(BColors.red(f"Error reading or parsing total bases file: {e}"), file=sys.stderr)
        sys.exit(1)

def load_and_filter_hits(input_files, pid_cutoff_fraction, pid_type):
    """
    Reads all input TSV files, filters by the chosen PID type, and returns a dictionary of
    all passing hits grouped by query_id.
    """
    hits_by_query = defaultdict(list)
    pid_cutoff_percent = pid_cutoff_fraction * 100.0
    pid_column = f"{pid_type}_pid"
    print(BColors.cyan(f"--- Processing {len(input_files)} input file(s) with {pid_column} >= {pid_cutoff_percent:.2f}%... ---"))
    
    for file_path in input_files:
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    try:
                        if float(row[pid_column]) >= pid_cutoff_percent:
                            hits_by_query[row['query_id']].append(row)
                    except (ValueError, KeyError):
                        continue
        except FileNotFoundError:
            print(BColors.yellow(f"Warning: Input file not found, skipping: {file_path}"), file=sys.stderr)
    return hits_by_query

def group_hits_by_final_category(hits_by_query, aro_to_family, consensus_threshold):
    """
    Resolves each read into its final category and groups the actual hit data
    by that category.
    """
    hits_by_final_category = defaultdict(list)
    for query_id, hits in hits_by_query.items():
        unique_aros_hit = {h['ARO_matched'] for h in hits}
        
        if len(unique_aros_hit) == 1:
            aro = list(unique_aros_hit)[0]
            family = aro_to_family.get(aro, "Unknown_Family")
            category_key = (family, aro)
        else:
            family_aro_pairs = [f"{aro_to_family.get(aro, 'Unknown_Family')};{aro}" for aro in unique_aros_hit]
            sorted_pairs_str = "|".join(sorted(family_aro_pairs))
            aro_description = f"multiple[{sorted_pairs_str}]"
            
            family_counts = Counter(aro_to_family.get(aro, "Unknown_Family") for aro in unique_aros_hit)
            most_common_family, count = family_counts.most_common(1)[0]
            family_description = most_common_family if count / len(unique_aros_hit) >= consensus_threshold else "multiple"
            category_key = (family_description, aro_description)
            
        hits_by_final_category[category_key].extend(hits)
        
    return hits_by_final_category

def calculate_coverage_for_hit_group(list_of_hits, aro_to_length):
    """
    Calculates the lateral coverage for a specific group of hits.
    For multiple-ARO groups, this returns the average coverage across all constituent AROs.
    """
    hits_by_aro = defaultdict(list)
    for hit in list_of_hits:
        hits_by_aro[hit['ARO_matched']].append(hit)

    coverages = []
    for aro, hits in hits_by_aro.items():
        coverage_set = set()
        for hit in hits:
            try:
                start = int(hit['position_on_ref'])
                length = int(hit['nucleotide_denominator'])
                coverage_set.update(range(start, start + length))
            except (ValueError, KeyError):
                continue
        
        total_len = aro_to_length.get(aro)
        if total_len and total_len > 0:
            coverages.append((len(coverage_set) / total_len) * 100.0)

    return sum(coverages) / len(coverages) if coverages else 0.0

def calculate_normalised_metrics(hits_by_final_category, aro_to_length, uscg_rpk, uscg_fpk, total_bases):
    """Calculates all metrics for each detailed AMR entry, including read and fragment counts."""
    print(BColors.cyan("--- Calculating metrics for detailed report... ---"))
    results = []
    total_bases_in_gb = (total_bases / 1e9) if total_bases > 0 else 0
    
    for category_key, list_of_hits in hits_by_final_category.items():
        family, aro = category_key
        
        all_read_ids = [h['query_id'] for h in list_of_hits]
        read_count = len(set(all_read_ids))
        
        fragment_ids = {re.sub(r'_\d+$', '', r_id) for r_id in all_read_ids}
        fragment_count = len(fragment_ids)
        
        coverage = calculate_coverage_for_hit_group(list_of_hits, aro_to_length)
        
        gene_length_bp = 0
        if aro.startswith("multiple[") or aro == "multiple":
            all_aros_in_group = {h['ARO_matched'] for h in list_of_hits}
            lengths = [aro_to_length.get(a, 0) for a in all_aros_in_group]
            valid_lengths = [l for l in lengths if l > 0]
            if valid_lengths:
                gene_length_bp = sum(valid_lengths) / len(valid_lengths)
        else:
            gene_length_bp = aro_to_length.get(aro, 0)

        rpk = (read_count / (gene_length_bp / 1000.0)) if gene_length_bp > 0 else 0.0
        rpkg = (rpk / total_bases_in_gb) if total_bases_in_gb > 0 else 0.0
        rpkpc_val, rpkpmc_val = "NA", "NA"
        if uscg_rpk is not None and uscg_rpk > 1e-9:
            rpkpc = rpk / uscg_rpk
            rpkpmc = rpkpc * 1_000_000
            rpkpc_val = f"{rpkpc:.4f}"
            rpkpmc_val = f"{rpkpmc:.2f}"
        elif uscg_rpk is not None and rpk == 0.0:
            rpkpc_val, rpkpmc_val = "0.0000", "0.00"

        fpk = (fragment_count / (gene_length_bp / 1000.0)) if gene_length_bp > 0 else 0.0
        fpkg = (fpk / total_bases_in_gb) if total_bases_in_gb > 0 else 0.0
        fpkpc_val, fpkpmc_val = "NA", "NA"

        if uscg_fpk is not None and uscg_fpk > 1e-9:
            fpkpc = fpk / uscg_fpk
            fpkpmc = fpkpc * 1_000_000
            fpkpc_val = f"{fpkpc:.4f}"
            fpkpmc_val = f"{fpkpmc:.2f}"
        elif uscg_fpk is not None and fpk == 0.0:
            fpkpc_val, fpkpmc_val = "0.0000", "0.00"

        results.append({
            "key": f"{family};{aro}",
            "Read_Count": read_count,
            "Fragment_Count": fragment_count,
            "Lateral_Coverage_%": f"{coverage:.2f}",
            "Gene_Length_bp": f"{gene_length_bp:.2f}" if isinstance(gene_length_bp, float) else str(gene_length_bp),
            "RPK": f"{rpk:.4f}", "RPKG": f"{rpkg:.4f}", "RPKPC": rpkpc_val, "RPKPMC": rpkpmc_val,
            "FPK": f"{fpk:.4f}", "FPKG": f"{fpkg:.4f}", "FPKPC": fpkpc_val, "FPKPMC": fpkpmc_val
        })
    results.sort(key=lambda x: x['key'])
    return results

def aggregate_for_clean_summary(hits_by_final_category):
    """Aggregates detailed hit groups into 'clean' groups."""
    clean_hit_groups = defaultdict(list)
    for (family, aro), list_of_hits in hits_by_final_category.items():
        if aro.startswith("multiple["):
            clean_key = (family, "multiple")
        else:
            clean_key = (family, aro)
        clean_hit_groups[clean_key].extend(list_of_hits)
    return clean_hit_groups

def write_summary_file(results, output_path):
    """Writes a normalised summary file (either detailed or clean)."""
    print(BColors.cyan(f"--- Writing normalised summary to: {output_path} ---"))
    try:
        with open(output_path, 'w', newline='') as f:
            header = [
                "#AMR_Gene_Family;ARO", "Read_Count", "Fragment_Count", "Lateral_Coverage_%",
                "Gene_Length_bp", "RPK", "FPK", "RPKG", "FPKG", "RPKPC", "FPKPC",
                "RPKPMC", "FPKPMC"
            ]
            writer = csv.DictWriter(f, fieldnames=header, delimiter='\t')
            
            f.write('\t'.join(header) + '\n')
            
            for row_data in results:
                row_to_write = {
                    "#AMR_Gene_Family;ARO": row_data["key"],
                    "Read_Count": row_data["Read_Count"],
                    "Fragment_Count": row_data["Fragment_Count"],
                    "Lateral_Coverage_%": row_data["Lateral_Coverage_%"],
                    "Gene_Length_bp": row_data["Gene_Length_bp"],
                    "RPK": row_data["RPK"], "RPKG": row_data["RPKG"],
                    "RPKPC": row_data["RPKPC"], "RPKPMC": row_data["RPKPMC"],
                    "FPK": row_data["FPK"], "FPKG": row_data["FPKG"],
                    "FPKPC": row_data["FPKPC"], "FPKPMC": row_data["FPKPMC"]
                }
                writer.writerow(row_to_write)
        print(BColors.green(f"--- Successfully wrote report: {output_path} ---"))
    except IOError as e:
        print(BColors.red(f"Error writing report file: {e}"), file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description="Combines, tabulates, and normalizes AMR homology hits, producing both a detailed and a clean summary."
    )
    parser.add_argument("-i", "--input-files", required=True, help="Comma-delimited list of input `_hits.tsv` files.")
    parser.add_argument("--metadata", required=True, help="Path to the AMR database metadata file.")
    parser.add_argument("--uscg-report", required=True, help="Path to the USCG quantification report.")
    parser.add_argument("--total-bases-file", required=True, help="Path to a file containing the total number of bases.")
    parser.add_argument("--tmp-dir", default=".", help="Directory to store the output report files.")
    parser.add_argument("--output-prefix", required=True, help="Prefix for the output summary files.")
    parser.add_argument("--pid-cutoff", type=float, default=0.9, help="Minimum percent identity to consider a hit (0.0-1.0 scale).")
    parser.add_argument("--pid-type", choices=['protein', 'nucleotide'], default='protein', help="PID type to use for filtering. Default: protein")
    parser.add_argument("--consensus", type=float, default=0.9, help="Minimum fraction of non-unique hits that must map to the same gene family to reach a consensus.")
    args = parser.parse_args()

    if not 0.0 <= args.pid_cutoff <= 1.0 or not 0.0 <= args.consensus <= 1.0:
        print(BColors.red("Error: --pid-cutoff and --consensus must be fractions between 0.0 and 1.0."), file=sys.stderr)
        sys.exit(1)
        
    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.tmp_dir, exist_ok=True)
    
    aro_to_family, aro_to_length = load_metadata(args.metadata)
    uscg_rpk, uscg_fpk = load_uscg_metrics(args.uscg_report)
    total_bases = load_total_bases(args.total_bases_file)
    
    hits_by_query = load_and_filter_hits(input_files_list, args.pid_cutoff, args.pid_type)
    print(BColors.green(f"--- Filtered and processed {len(hits_by_query)} unique reads from all files. ---"))

    hits_by_final_category = group_hits_by_final_category(hits_by_query, aro_to_family, args.consensus)
    print(BColors.green(f"--- Aggregated reads into {len(hits_by_final_category)} detailed categories. ---"))

    detailed_normalised_results = calculate_normalised_metrics(hits_by_final_category, aro_to_length, uscg_rpk, uscg_fpk, total_bases)
    detailed_output_path = os.path.join(args.tmp_dir, f"{args.output_prefix}_homscan_detailed.tsv")
    write_summary_file(detailed_normalised_results, detailed_output_path)

    clean_hit_groups = aggregate_for_clean_summary(hits_by_final_category)
    print(BColors.green(f"--- Aggregated detailed results into {len(clean_hit_groups)} clean categories. ---"))

    clean_normalised_results = calculate_normalised_metrics(clean_hit_groups, aro_to_length, uscg_rpk, uscg_fpk, total_bases)
    clean_output_path = os.path.join(args.tmp_dir, f"{args.output_prefix}_homscan.tsv")
    write_summary_file(clean_normalised_results, clean_output_path)
    
    print(BColors.green("\n\n--- All Processing Complete ---"))
    print("The following output files were generated:")
    print(f"- {detailed_output_path}")
    print(f"- {clean_output_path}")

if __name__ == "__main__":
    main()