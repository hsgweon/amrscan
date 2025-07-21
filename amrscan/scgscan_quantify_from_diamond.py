#!/usr/bin/env python3

import argparse
import os
import sys
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

def load_gene_lengths(filepath):
    """Loads a two-column TSV file (gene_id, length) into a dictionary."""
    print(BColors.cyan(f"--- Loading gene lengths from: {filepath} ---"))
    if not os.path.exists(filepath):
        print(BColors.red(f"Error: Gene lengths file not found at '{filepath}'"), file=sys.stderr)
        return None
    
    lengths_map = {}
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline()
            if not (first_line.startswith("USCG_ID") or first_line.startswith("#")):
                f.seek(0)

            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    lengths_map[parts[0]] = int(parts[1])
    except (IOError, ValueError) as e:
        print(BColors.red(f"Error reading or parsing gene lengths file: {e}"), file=sys.stderr)
        return None
    
    print(BColors.green(f"--- Loaded lengths for {len(lengths_map)} genes."))
    return lengths_map

def process_diamond_file(diamond_path, evalue_cutoff):
    """
    Processes a DIAMOND tabular output file to count hits per reference gene.
    """
    print(BColors.cyan(f"--- Processing DIAMOND file: {diamond_path} ---"))
    counts = Counter()
    total_alignments = 0
    passing_alignments = 0
    
    try:
        with open(diamond_path, 'r') as f:
            for line in f:
                total_alignments += 1
                fields = line.strip().split('\t')
                
                if len(fields) >= 12:
                    full_gene_id = fields[1]
                    evalue = float(fields[10])
                    
                    if evalue <= evalue_cutoff:
                        # <<< THIS IS THE FIX: Extract the base gene name >>>
                        # Splits 'rpoB___409' into 'rpoB' and '409', takes the first part.
                        base_gene_id = full_gene_id.split('___')[0]
                        
                        counts[base_gene_id] += 1
                        passing_alignments += 1

    except FileNotFoundError:
        print(BColors.red(f"Error: DIAMOND file not found at '{diamond_path}'"), file=sys.stderr)
        return None, 0, 0
    
    print(BColors.green(f"--- Processed {total_alignments} alignments. Found {passing_alignments} hits passing E-value cutoff for {len(counts)} unique genes."))
    return counts, total_alignments, passing_alignments

def write_report(output_path, uscg_counts, uscg_lengths, total_alignments):
    """Calculates metrics and writes the final summary report."""
    print(BColors.cyan(f"--- Calculating metrics and writing report to: {output_path} ---"))
    
    total_mapped_frags = sum(uscg_counts.values())
    total_ref_len = sum(uscg_lengths.values())
    uscg_rpks = {}

    for gene, length in uscg_lengths.items():
        count = uscg_counts.get(gene, 0)
        uscg_rpks[gene] = (count / (length / 1000.0)) if length > 0 else 0.0

    overall_rpk = (total_mapped_frags / (total_ref_len / 1000.0)) if total_ref_len > 0 else 0.0

    try:
        with open(output_path, 'w') as f:
            f.write("### USCG Quantification Summary ###\n")
            f.write("USCG_ID\tMapped_Fragment_Count\tRPK\n")
            
            for gene_id in sorted(uscg_lengths.keys()):
                count = uscg_counts.get(gene_id, 0)
                rpk = uscg_rpks.get(gene_id, 0.0)
                f.write(f"{gene_id}\t{count}\t{rpk:.4f}\n")

            f.write("\n### Overall Sample Statistics ###\n")
            f.write(f"Total_Alignments_In_DIAMOND_File\t{total_alignments}\n")
            f.write(f"Total_Fragments_Mapped_To_Any_USCG_Passing_Filters\t{total_mapped_frags}\n")
            f.write(f"Total_Reference_USCG_Length_AA_Combined\t{total_ref_len}\n")
            f.write(f"Overall_RPK_Across_All_USCGs\t{overall_rpk:.4f}\n")
        
        print(BColors.green("--- Successfully wrote report."))
    except IOError as e:
        print(BColors.red(f"Error writing report file: {e}"), file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description="Quantifies Universal Single-Copy Genes (USCGs) from a DIAMOND tabular output file."
    )
    parser.add_argument(
        "-i", "--diamond-files",
        required=True,
        help="Comma-delimited list of input DIAMOND tabular files."
    )
    parser.add_argument(
        "-l", "--gene-lengths",
        required=True,
        help="Path to a tab-separated file containing USCG IDs and their lengths in AMINO ACIDS (USCG_ID\\tLength)."
    )
    parser.add_argument(
        "-o", "--output-report",
        required=True,
        help="Path for the output summary TSV report."
    )
    parser.add_argument(
        "-e", "--evalue",
        type=float,
        default=1e-5,
        help="E-value cutoff for considering an alignment a valid hit. (default: 1e-5)"
    )
    args = parser.parse_args()

    uscg_lengths = load_gene_lengths(args.gene_lengths)
    if uscg_lengths is None:
        sys.exit(1)

    diamond_files_list = [f.strip() for f in args.diamond_files.split(',') if f.strip()]
    if not diamond_files_list:
        print(BColors.red("Error: No input DIAMOND files provided."), file=sys.stderr)
        sys.exit(1)

    aggregated_counts = Counter()
    total_alignments_processed = 0

    for diamond_file in diamond_files_list:
        counts, alignments_in_file, _ = process_diamond_file(diamond_file, args.evalue)
        if counts is None:
            continue
        
        aggregated_counts.update(counts)
        total_alignments_processed += alignments_in_file

    write_report(args.output_report, aggregated_counts, uscg_lengths, total_alignments_processed)
    
    print(BColors.green("\n\n--- USCG Quantification from DIAMOND Complete ---"))

if __name__ == "__main__":
    main()