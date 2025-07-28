#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
import shutil

class BColors:
    """A helper class to add color to terminal output."""
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


def check_dependencies():
    """Checks if DIAMOND is in the system's PATH."""
    print(BColors.cyan("--- Checking for required dependencies (diamond)..."))
    
    if not shutil.which("diamond"):
        print(BColors.red("Error: 'diamond' is not found in your PATH. Please install it."), file=sys.stderr)
        sys.exit(1)
        
    print(BColors.green("--- Dependencies found."))

def run_command(command, description):
    """Executes a simple command and checks for errors."""
    print(BColors.cyan(f"\n--- {description} ---"))
    cmd_str = ' '.join(command)
    print(f"Executing: {BColors.yellow(cmd_str)}")
    
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        if result.stderr:
            print(BColors.cyan("--- DIAMOND run statistics (from stderr) ---"))
            print(result.stderr.strip())
    except subprocess.CalledProcessError as e:
        print(BColors.red(f"\nError executing command: {cmd_str}"), file=sys.stderr)
        print(BColors.red(f"Return code: {e.returncode}"), file=sys.stderr)
        print(BColors.red("\n--- STDOUT ---"), file=sys.stderr)
        print(e.stdout, file=sys.stderr)
        print(BColors.red("\n--- STDERR ---"), file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="A wrapper to map nucleotide reads to a PROTEIN reference database using DIAMOND, producing a tabular output of the single best hit for each input file."
    )
    parser.add_argument(
        "-i", "--input-files",
        required=True,
        help="Comma-delimited list of input FASTQ/FASTA files (each file is processed independently)."
    )
    parser.add_argument(
        "-d", "--db",
        required=True,
        help="Path to the reference PROTEIN database FASTA file."
    )
    parser.add_argument(
        "--tmp-dir",
        default=".",
        help="Directory for storing the DIAMOND index and output files. Defaults to the current directory."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=8,
        help="Number of threads to use for the alignment step. (default: 8)"
    )
    parser.add_argument(
        "--sensitivity",
        choices=['fast', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive'],
        default='fast',
        help="DIAMOND sensitivity mode. (default: fast)"
    )
    
    args = parser.parse_args()

    check_dependencies()
    
    os.makedirs(args.tmp_dir, exist_ok=True)
    
    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    db_basename = os.path.basename(args.db)
    db_index_path = os.path.join(args.tmp_dir, f"{db_basename}.dmnd")
    
    # Build DIAMOND database index
    diamond_db_cmd = ['diamond', 'makedb', '--in', args.db, '--db', db_index_path, '--threads', str(args.threads)]
    run_command(diamond_db_cmd, "Building DIAMOND database index")

    # Alignment Step: Loop through each input file
    created_tsv_files = []
    for input_file in input_files_list:
        input_basename = os.path.basename(input_file)
        output_tsv_path = os.path.join(args.tmp_dir, f"{input_basename}.diamond.tsv")

        aligner_cmd = [
            'diamond', 'blastx',
            '--db', db_index_path,
            '--threads', str(args.threads),
            f'--{args.sensitivity}',
            '--query', input_file,
            '--out', output_tsv_path,
            '--max-target-seqs', '1',
            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
        ]
        
        run_command(aligner_cmd, f"Mapping reads for {input_basename} with DIAMOND")
        print(BColors.green(f"--- Successfully created DIAMOND output file: {output_tsv_path} ---"))
        created_tsv_files.append(output_tsv_path)

    print(BColors.green("\n\n--- All Processing Complete ---"))
    print("The following output files were generated:")
    for f in created_tsv_files:
        print(f"- {f}")

if __name__ == "__main__":
    main()