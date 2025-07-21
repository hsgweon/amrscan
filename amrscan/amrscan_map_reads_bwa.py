#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys
import shutil

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


def check_dependencies():
    """Checks if required command-line tools are in the system's PATH."""
    print(BColors.cyan("--- Checking for required dependencies (bwa, samtools)..."))
    
    if not shutil.which("bwa"):
        print(BColors.red("Error: 'bwa' is not found in your PATH. Please install it."), file=sys.stderr)
        sys.exit(1)
        
    if not shutil.which("samtools"):
        print(BColors.red("Error: 'samtools' is not found in your PATH. Please install it."), file=sys.stderr)
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
            print(result.stderr.strip())
    except subprocess.CalledProcessError as e:
        print(BColors.red(f"\nError executing command: {cmd_str}"), file=sys.stderr)
        print(BColors.red(f"Return code: {e.returncode}"), file=sys.stderr)
        print(BColors.red("\n--- STDOUT ---"), file=sys.stderr)
        print(e.stdout, file=sys.stderr)
        print(BColors.red("\n--- STDERR ---"), file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        sys.exit(1)

def run_pipeline(aligner_cmd, samtools_cmd, input_file):
    """Executes a shell pipeline (bwa | samtools) for a single input file."""
    print(BColors.cyan(f"\n--- Mapping reads for: {os.path.basename(input_file)} ---"))
    pipeline_str = f"{' '.join(aligner_cmd)} | {' '.join(samtools_cmd)}"
    print(f"Executing pipeline: {BColors.yellow(pipeline_str)}")

    try:
        aligner_process = subprocess.Popen(aligner_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        samtools_process = subprocess.Popen(samtools_cmd, stdin=aligner_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        aligner_process.stdout.close()
        
        samtools_stdout, samtools_stderr = samtools_process.communicate()
        aligner_stdout, aligner_stderr = aligner_process.communicate()

        if aligner_process.returncode != 0:
            print(BColors.red(f"\nError in BWA part of the pipeline for file {input_file}."), file=sys.stderr)
            print(BColors.red(f"BWA return code: {aligner_process.returncode}"), file=sys.stderr)
            print(BColors.red("\n--- BWA STDERR ---"), file=sys.stderr)
            print(aligner_stderr.decode('utf-8'), file=sys.stderr)
            sys.exit(1)
            
        if samtools_process.returncode != 0:
            print(BColors.red(f"\nError in samtools part of the pipeline for file {input_file}."), file=sys.stderr)
            print(BColors.red(f"samtools return code: {samtools_process.returncode}"), file=sys.stderr)
            print(BColors.red("\n--- SAMTOOLS STDERR ---"), file=sys.stderr)
            print(samtools_stderr.decode('utf-8'), file=sys.stderr)
            sys.exit(1)

        if aligner_stderr:
             print(BColors.cyan("--- BWA run statistics (from stderr) ---"))
             print(aligner_stderr.decode('utf-8').strip())

    except Exception as e:
        print(BColors.red(f"An unexpected error occurred while running the pipeline for {input_file}: {e}"), file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="A wrapper to map reads to a reference database using BWA."
    )
    parser.add_argument(
        "-i", "--input-files",
        required=True,
        help="Comma-delimited list of input FASTQ/FASTA files."
    )
    parser.add_argument(
        "-d", "--db",
        required=True,
        help="Path to the reference database FASTA file."
    )
    parser.add_argument(
        "--tmp-dir",
        default=".",
        help="Directory for storing temporary files (index, output SAMs). Defaults to the current directory."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=8,
        help="Number of threads to use for the alignment step. (default: 8)"
    )
    parser.add_argument(
        "--best-hit-only",
        action="store_true",
        help="If specified, BWA will only report the single best hit for each read (secondary alignments are discarded)."
    )
    
    args = parser.parse_args()

    check_dependencies()
    
    os.makedirs(args.tmp_dir, exist_ok=True)
    
    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    db_basename = os.path.basename(args.db)
    db_index_path = os.path.join(args.tmp_dir, db_basename)
    
    bwa_index_cmd = ['bwa', 'index', '-p', db_index_path, args.db]
    run_command(bwa_index_cmd, "Building BWA database index")

    created_files = []
    for input_file in input_files_list:
        input_basename = os.path.basename(input_file)
        output_sam_path = os.path.join(args.tmp_dir, f"{input_basename}.sam")

        if args.best_hit_only:
            print(BColors.yellow("\n--- Running in 'best-hit-only' mode. Secondary alignments will be discarded. ---"))
            # BWA command without '-a' to produce primary/secondary alignments
            aligner_cmd = ['bwa', 'mem', '-t', str(args.threads), db_index_path, input_file]
            # Samtools command to filter unmapped (4) AND secondary (256) alignments. 4+256=260
            samtools_cmd = ['samtools', 'view', '-F', '260', '-', '-o', output_sam_path]
        else:
            aligner_cmd = ['bwa', 'mem', '-t', str(args.threads), '-a', db_index_path, input_file]
            # Samtools command to only filter unmapped reads (4)
            samtools_cmd = ['samtools', 'view', '-F', '4', '-', '-o', output_sam_path]
        
        run_pipeline(aligner_cmd, samtools_cmd, input_file)
        print(BColors.green(f"--- Successfully created mapped reads file: {output_sam_path} ---"))
        created_files.append(output_sam_path)

    print(BColors.green("\n\n--- All Processing Complete ---"))
    print("The following output files were generated:")
    for f in created_files:
        print(f"- {f}")

if __name__ == "__main__":
    main()