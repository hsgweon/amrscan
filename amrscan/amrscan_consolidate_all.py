#!/usr/bin/env python3

import argparse
import os
import sys
import shutil
from pathlib import Path

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

def main():
    parser = argparse.ArgumentParser(
        description="Consolidates the final summary reports from the homscan and varscan pipelines into the main output directory."
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

    print(BColors.cyan(f"--- Consolidating final reports for prefix: {filename_prefix} ---"))

    files_to_copy = [
        ('homscan', f"{filename_prefix}_homscan.tsv"),
        ('homscan', f"{filename_prefix}_homscan_detailed.tsv"),
        ('varscan', f"{filename_prefix}_varscan.tsv")
    ]
    
    final_reports = []

    for sub_dir, filename in files_to_copy:
        source_path = main_output_dir / 'tmp' / sub_dir / filename
        destination_path = main_output_dir / filename

        if source_path.exists():
            try:
                shutil.copy(source_path, destination_path)
                print(BColors.green(f"Successfully copied: {filename}"))
                final_reports.append(str(destination_path))
            except IOError as e:
                print(BColors.red(f"Error copying {filename}: {e}"), file=sys.stderr)
        else:
            print(BColors.yellow(f"Warning: Source file not found, skipping: {source_path}"), file=sys.stderr)

            try:
                with open(destination_path, 'w') as f:
                    if "homscan" in filename:
                        f.write("#AMR_Gene_Family;ARO\tCount\tLateral_Coverage_%\tGene_Length_bp\tRPK\tRPKG\tRPKPC\tRPKPMC\n")
                    elif "varscan" in filename:
                        f.write("#AMR_Gene_Family;ARO\tCount\tGene_Length_bp\tRPK\tRPKG\tRPKPC\tRPKPMC\n")

                print(BColors.yellow(f"Created an empty placeholder for: {filename}"), file=sys.stderr)
                final_reports.append(str(destination_path))
            except IOError as e:
                print(BColors.red(f"Error creating placeholder file {destination_path}: {e}"), file=sys.stderr)


    print(BColors.green("\n\n--- Consolidation Complete ---"))
    print("The following final summary reports are available in the main output directory:")
    for report in sorted(final_reports):
        print(f"- {report}")

if __name__ == "__main__":
    main()