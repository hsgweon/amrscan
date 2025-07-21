#!/usr/bin/env python3

import argparse
import os
import sys
import csv
from collections import Counter

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

    aro_to_family = {}
    aro_to_length = {}
    try:
        with open(metadata_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            if 'SeqNucLength' not in reader.fieldnames:
                print(BColors.red("Error: Metadata file must contain a 'SeqNucLength' column for normalisation."), file=sys.stderr)
                sys.exit(1)
            for row in reader:
                aro_to_family[row['Sequence_ID']] = row['AMR_Gene_Family']
                try:
                    aro_to_length[row['Sequence_ID']] = int(row['SeqNucLength'])
                except (ValueError, TypeError):
                    aro_to_length[row['Sequence_ID']] = 0
    except (Exception, KeyError) as e:
        print(BColors.red(f"Error reading or parsing metadata file: {e}"), file=sys.stderr)
        sys.exit(1)

    print(BColors.green(f"--- Loaded metadata for {len(aro_to_family)} AROs."))
    return aro_to_family, aro_to_length

def load_uscg_rpk(uscg_report_path):
    """Loads the Overall_RPK_Across_All_USCGs value from the USCG report."""
    print(BColors.cyan(f"--- Loading USCG RPK from: {uscg_report_path} ---"))
    if not os.path.exists(uscg_report_path):
        print(BColors.red(f"Error: USCG report file not found at '{uscg_report_path}'"), file=sys.stderr)
        return None
    try:
        with open(uscg_report_path, 'r') as f:
            for line in f:
                if line.startswith("Overall_RPK_Across_All_USCGs"):
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        rpk_val = float(parts[1])
                        print(BColors.green(f"--- Found Overall USCG RPK: {rpk_val:.4f} ---"))
                        return rpk_val
    except (IOError, ValueError) as e:
        print(BColors.red(f"Error reading or parsing USCG report file: {e}"), file=sys.stderr)
        return None
    print(BColors.red("Error: 'Overall_RPK_Across_All_USCGs' not found in the USCG report."), file=sys.stderr)
    return None

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

def count_variant_hits(input_files, pid_cutoff_fraction):
    """
    Reads all input _variant_hits.tsv files, filters by PID, and counts confirmed
    reads for each variant ARO.
    """
    final_counts = Counter()
    pid_cutoff_percent = pid_cutoff_fraction * 100.0
    
    print(BColors.cyan(f"--- Processing {len(input_files)} variant hits file(s) with PID >= {pid_cutoff_percent:.2f}%... ---"))
    for file_path in input_files:
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    try:
                        if float(row['nucleotide_pid']) >= pid_cutoff_percent:
                            final_counts[row['ARO_matched']] += 1
                    except (ValueError, KeyError):
                        print(BColors.yellow(f"Warning: Skipping row in {file_path} due to missing or invalid columns."), file=sys.stderr)
                        continue
        except FileNotFoundError:
            print(BColors.yellow(f"Warning: Input file not found, skipping: {file_path}"), file=sys.stderr)
    
    return final_counts

def calculate_and_write_summary(final_counts, aro_to_family, aro_to_length, uscg_rpk, total_bases, output_path):
    """Calculates all metrics and writes the final normalised variant summary."""
    print(BColors.cyan(f"--- Writing normalised variant summary to: {output_path} ---"))
    total_bases_in_gb = (total_bases / 1e9) if total_bases > 0 else 0
    
    results = []
    for aro, count in final_counts.items():
        family = aro_to_family.get(aro, "Unknown_Family")
        gene_length_bp = aro_to_length.get(aro, 0)

        rpk = (count / (gene_length_bp / 1000.0)) if gene_length_bp > 0 else 0.0
        rpkg = (rpk / total_bases_in_gb) if total_bases_in_gb > 0 else 0.0
        
        rpkpc_val, rpkpmc_val = "NA", "NA"
        if uscg_rpk is not None and uscg_rpk > 1e-9:
            rpkpc = rpk / uscg_rpk
            rpkpmc = rpkpc * 1_000_000
            rpkpc_val = f"{rpkpc:.4f}"
            rpkpmc_val = f"{rpkpmc:.2f}"
        elif uscg_rpk is not None and rpk == 0.0:
            rpkpc_val, rpkpmc_val = "0.0000", "0.00"

        results.append({
            "key": f"{family};{aro}", "Count": count,
            "Gene_Length_bp": str(gene_length_bp),
            "RPK": f"{rpk:.4f}", "RPKG": f"{rpkg:.4f}", "RPKPC": rpkpc_val, "RPKPMC": rpkpmc_val
        })
    
    results.sort(key=lambda x: x['key'])

    try:
        with open(output_path, 'w', newline='') as f:
            header = ["#AMR_Gene_Family;ARO", "Count", "Gene_Length_bp", "RPK", "RPKG", "RPKPC", "RPKPMC"]
            writer = csv.DictWriter(f, fieldnames=header, delimiter='\t')
            f.write('\t'.join(header) + '\n')
            for row_data in results:
                row_to_write = {
                    "#AMR_Gene_Family;ARO": row_data["key"], "Count": row_data["Count"],
                    "Gene_Length_bp": row_data["Gene_Length_bp"], "RPK": row_data["RPK"],
                    "RPKG": row_data["RPKG"], "RPKPC": row_data["RPKPC"], "RPKPMC": row_data["RPKPMC"]
                }
                writer.writerow(row_to_write)
        print(BColors.green(f"--- Successfully wrote report: {output_path} ---"))
    except IOError as e:
        print(BColors.red(f"Error writing report file: {e}"), file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description="Tabulates and normalises confirmed variant hits into a summary report."
    )
    parser.add_argument("-i", "--input-files", required=True, help="Comma-delimited list of input `_variant_hits.tsv` files.")
    parser.add_argument("--metadata", required=True, help="Path to the AMR database metadata file (TSV format, must include 'SeqNucLength' column).")
    parser.add_argument("--uscg-report", required=True, help="Path to the USCG quantification report containing the overall RPK.")
    parser.add_argument("--total-bases-file", required=True, help="Path to a file containing the total number of bases for RPKG calculation.")
    parser.add_argument("--tmp-dir", default=".", help="Directory to store the output report file. Defaults to the current directory.")
    parser.add_argument("--output-prefix", required=True, help="Prefix for the output summary file (e.g., 'MyProject').")
    parser.add_argument(
        "--pid-cutoff",
        type=float,
        default=0.95,
        help="Minimum nucleotide PID to consider a hit (0.0-1.0 scale). Default: 0.95"
    )
    
    args = parser.parse_args()
        
    if not 0.0 <= args.pid_cutoff <= 1.0:
        print(BColors.red("Error: --pid-cutoff must be a fraction between 0.0 and 1.0."), file=sys.stderr)
        sys.exit(1)

    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.tmp_dir, exist_ok=True)
    
    aro_to_family, aro_to_length = load_metadata(args.metadata)
    uscg_rpk = load_uscg_rpk(args.uscg_report)
    total_bases = load_total_bases(args.total_bases_file)
    
    final_counts = count_variant_hits(input_files_list, args.pid_cutoff)
    print(BColors.green(f"--- Aggregated results for {len(final_counts)} unique variants passing filters. ---"))

    summary_output_path = os.path.join(args.tmp_dir, f"{args.output_prefix}_varscan.tsv")
    calculate_and_write_summary(final_counts, aro_to_family, aro_to_length, uscg_rpk, total_bases, summary_output_path)
    
    print(BColors.green("\n\n--- Variant Tabulation and Normalization Complete ---"))

if __name__ == "__main__":
    main()