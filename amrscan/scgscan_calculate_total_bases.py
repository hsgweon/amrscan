#!/usr/bin/env python3

# amrscan_pipeline/amrscan/scgscan_calculate_total_bases.py
import argparse
import os
import sys
import gzip
import bz2

class BColors:
    """A helper class to add color to terminal output."""
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
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

def open_fastq_file(filepath):
    """Opens a FASTQ file, transparently handling .gz, .bz2, or uncompressed."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    elif filepath.endswith('.bz2'):
        return bz2.open(filepath, 'rt')
    else:
        return open(filepath, 'r')

def main():
    parser = argparse.ArgumentParser(
        description="Calculates the total number of bases from the sequence lines of FASTQ files."
    )
    parser.add_argument(
        "-i", "--input-files",
        required=True,
        help="Comma-delimited list of input FASTQ files (can be .gz, .bz2, or uncompressed)."
    )
    parser.add_argument(
        "-o", "--output-file",
        required=True,
        help="Path to the output file where the total base count will be written."
    )
    args = parser.parse_args()

    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    total_bases = 0
    for filepath in input_files_list:
        if not os.path.exists(filepath):
            print(BColors.red(f"Error: Input file not found: {filepath}"), file=sys.stderr)
            sys.exit(1)
        
        print(BColors.cyan(f"--- Processing file: {os.path.basename(filepath)} ---"))
        try:
            with open_fastq_file(filepath) as f:
                line_num = 0
                for line in f:
                    line_num += 1
                    if line_num % 4 == 2:
                        total_bases += len(line.strip())
        except Exception as e:
            print(BColors.red(f"Error processing file {filepath}: {e}"), file=sys.stderr)
            sys.exit(1)

    print(BColors.cyan(f"--- Total bases calculated: {total_bases:,} ---"))

    try:
        output_dir = os.path.dirname(args.output_file)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            
        with open(args.output_file, 'w') as f_out:
            f_out.write(str(total_bases))
        print(BColors.green(f"--- Successfully wrote total base count to: {args.output_file} ---"))
    except IOError as e:
        print(BColors.red(f"Error writing to output file {args.output_file}: {e}"), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()