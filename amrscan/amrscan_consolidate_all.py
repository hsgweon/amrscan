#!/usr/bin/env python3

# amrscan_pipeline/amrscan/amrscan_consolidate_all.py
import argparse
import os
import sys
import re
from pathlib import Path
import pandas as pd

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

FINAL_HEADERS = {
    "homscan": [
        "AMR_Gene_Family", "ARO", "CARD_Short_Name_Misc", "Read_Count", "Fragment_Count", 
        "Lateral_Coverage_%", "Gene_Length_bp", "RPK", "FPK", "RPKG", "FPKG", 
        "RPKPC", "FPKPC", "RPKPMC", "FPKPMC", "AMR_Gene_Family;ARO"
    ],
    "varscan": [
        "AMR_Gene_Family", "ARO", "CARD_Short_Name_Misc", "Read_Count", "Fragment_Count", 
        "Gene_Length_bp", "RPK", "FPK", "RPKG", "FPKG", "RPKPC", "FPKPC", "RPKPMC", "FPKPMC",
        "AMR_Gene_Family;ARO"
    ],
    "homscan_MAP": [
        "AMR_Gene_Family", "ARO", "CARD_Short_Name_Misc", "Read_Count", "Fragment_Count", 
        "Lateral_Coverage_%", "Gene_Length_bp", "RPK", "FPK", "RPKG", "FPKG", 
        "RPKPC", "FPKPC", "RPKPMC", "FPKPMC", "Resolved_RPKG", "Top_ARO", "Allocation_Proportions",
        "AMR_Gene_Family;ARO"
    ]
}

def get_report_type(filename):
    """Determines the report type from the filename."""
    if "homscan_MAP.tsv" in filename:
        return "homscan_MAP"
    if "homscan.tsv" in filename:
        return "homscan"
    if "varscan.tsv" in filename:
        return "varscan"
    return None

def parse_amr_key(key_string):
    """
    Parses the '#AMR_Gene_Family;ARO' key into three parts.
    Handles complex family names containing semicolons by splitting on the last
    semicolon that precedes an ARO identifier.
    """
    match = re.search(r';(ARO_\d+.*)', key_string)

    if match:
        split_point = match.start()
        family = key_string[:split_point]
        aro_part = match.group(1)
    else:
        try:
            family, aro_part = key_string.split(';', 1)
        except ValueError:
            return key_string, "N/A", "N/A"

    if aro_part.startswith("multiple"):
        return family, "multiple", "N/A"
    
    if '__' in aro_part:
        aro_id_block, misc_block = aro_part.split('__', 1)
        return family, aro_id_block, misc_block
    else:
        return family, aro_part, "N/A"

def reformat_and_write_summary(source_path, destination_path, report_type):
    """
    Reads a summary file, reformats it to be user-friendly, and writes the result.
    Creates an empty placeholder with a hashed header if the source does not exist.
    """
    header = FINAL_HEADERS.get(report_type)
    if not header:
        print(BColors.red(f"Error: Unknown report type for {source_path.name}. Cannot reformat."), file=sys.stderr)
        return

    if not source_path.exists() or os.path.getsize(source_path) == 0:
        print(BColors.yellow(f"Warning: Source file not found or empty, creating placeholder: {source_path}"), file=sys.stderr)
        try:
            with open(destination_path, 'w') as f:
                f.write('#' + '\t'.join(header) + '\n')
            print(BColors.yellow(f"Created an empty placeholder for: {destination_path.name}"), file=sys.stderr)
        except IOError as e:
            print(BColors.red(f"Error creating placeholder file {destination_path}: {e}"), file=sys.stderr)
        return

    try:
        df = pd.read_csv(source_path, sep='\t')
        
        if df.empty:
            with open(destination_path, 'w') as f:
                f.write('#' + '\t'.join(header) + '\n')
            return

        key_col = '#AMR_Gene_Family;ARO'
        if key_col in df.columns:
            df[['AMR_Gene_Family', 'ARO', 'CARD_Short_Name_Misc']] = df[key_col].apply(
                lambda x: pd.Series(parse_amr_key(x))
            )
        else:
             print(BColors.red(f"Error: Key column '{key_col}' not found in {source_path.name}."), file=sys.stderr)
             return

        final_df = pd.DataFrame()
        for col in header:
            if col in df.columns:
                final_df[col] = df[col]
            else:
                final_df[col] = None
        
        with open(destination_path, 'w', newline='') as f_out:
            f_out.write('#' + '\t'.join(final_df.columns) + '\n')
            final_df.to_csv(f_out, sep='\t', index=False, header=False, float_format='%.6f', na_rep='N/A')

        print(BColors.green(f"Successfully reformatted and copied: {destination_path.name}"))

    except Exception as e:
        print(BColors.red(f"Error processing file {source_path.name}: {e}"), file=sys.stderr)
        with open(destination_path, 'w') as f:
            f.write('#' + '\t'.join(header) + '\n')

def main():
    parser = argparse.ArgumentParser(
        description="Consolidates and reformats the final summary reports from the homscan and varscan pipelines."
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="The common prefix used for all output files and directories (e.g., 'SRR10842871'). This should be the full path to the output directory."
    )
    args = parser.parse_args()

    main_output_dir = Path(args.output_prefix)
    filename_prefix = main_output_dir.name
    
    if not main_output_dir.is_dir():
        print(BColors.red(f"Error: Main output directory '{main_output_dir}' does not exist. Please ensure previous steps ran correctly."), file=sys.stderr)
        sys.exit(1)

    print(BColors.cyan(f"--- Consolidating and reformatting final reports for prefix: {filename_prefix} ---"))

    files_to_process = [
        ('homscan', f"{filename_prefix}_homscan.tsv"),
        ('homscan', f"{filename_prefix}_homscan_MAP.tsv"),
        ('varscan', f"{filename_prefix}_varscan.tsv")
    ]
    
    final_reports = []

    for sub_dir, filename in files_to_process:
        source_path = main_output_dir / 'tmp' / sub_dir / filename
        destination_path = main_output_dir / filename
        report_type = get_report_type(filename)
        
        if report_type:
            reformat_and_write_summary(source_path, destination_path, report_type)
            final_reports.append(str(destination_path))
        else:
            print(BColors.yellow(f"Warning: Could not determine report type for '{filename}'. Skipping."), file=sys.stderr)

    print(BColors.green("\n\n--- Consolidation Complete ---"))
    print("The following final summary reports are available in the main output directory:")
    for report in sorted(final_reports):
        print(f"- {report}")

if __name__ == "__main__":
    main()