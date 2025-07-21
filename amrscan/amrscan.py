#!/usr/bin/env python3

import argparse
import logging
import subprocess
import sys
import time
from pathlib import Path
import os
import platform
import re
import shutil
import glob

try:
    import importlib.metadata
    from importlib import resources
except ImportError:
    import pkg_resources

__version__ = "1.0.0"

# ANSI escape codes for colors
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    GREY = '\033[90m'

class ColoredFormatter(logging.Formatter):
    """Custom formatter to add colors to log messages for console output."""
    LOG_FORMAT = "%(message)s"
    FORMATS = {
        logging.DEBUG: Colors.GREY + LOG_FORMAT + Colors.ENDC,
        logging.INFO: Colors.GREEN + LOG_FORMAT + Colors.ENDC,
        logging.WARNING: Colors.WARNING + LOG_FORMAT + Colors.ENDC,
        logging.ERROR: Colors.FAIL + Colors.BOLD + LOG_FORMAT + Colors.ENDC,
        logging.CRITICAL: Colors.FAIL + Colors.BOLD + Colors.UNDERLINE + LOG_FORMAT + Colors.ENDC,
    }
    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

def get_py_dependency_versions():
    """Retrieves the versions of specified Python packages."""
    dependencies = ['pip', 'setuptools']
    versions = {}
    use_importlib = 'importlib' in sys.modules and hasattr(sys.modules.get('importlib'), 'metadata')

    for package in dependencies:
        try:
            if use_importlib:
                versions[package] = importlib.metadata.version(package)
            else:
                versions[package] = pkg_resources.get_distribution(package).version
        except Exception:
            versions[package] = "Not Found"
    return versions

def get_tool_versions():
    """Retrieves the versions of external command-line tools."""
    versions = {}
    try:
        result = subprocess.run(['bwa'], capture_output=True, text=True, check=False)
        version_line = next((line for line in result.stderr.splitlines() if line.startswith('Version:')), None)
        versions['BWA'] = version_line.split(':', 1)[1].strip() if version_line else 'Unknown'
    except Exception:
        versions['BWA'] = 'Not Found'
    try:
        result = subprocess.run(['diamond', '--version'], capture_output=True, text=True, check=True)
        versions['Diamond'] = result.stdout.strip().split()[-1]
    except Exception:
        versions['Diamond'] = 'Not Found'
    return versions

def setup_logging(log_file):
    """Sets up logging to a file (plain) and console (colored)."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers = []
    file_handler = logging.FileHandler(log_file, mode='w')
    file_formatter = logging.Formatter("%(asctime)s [%(levelname)s] - %(message)s")
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(ColoredFormatter())
    logger.addHandler(console_handler)

def log_boxed_header(title, color_class):
    """Logs a formatted, colored, boxed section header."""
    padding = 2
    width = len(title) + padding * 2
    border_line = f"+{'-' * width}+"
    title_line = f"|{' ' * padding}{title.upper()}{' ' * padding}|"
    
    header_text = f"{border_line}\n{title_line}\n{border_line}"
    
    print(f"\n{color_class}{Colors.BOLD}{header_text}{Colors.ENDC}\n")
    
    logging.getLogger().handlers[0].handle(logging.makeLogRecord({
        'msg': f"\n--- {title.upper()} ---\n", 'levelname': 'INFO', 'levelno': logging.INFO,
        'pathname': '', 'lineno': 0, 'exc_info': None, 'func': None
    }))

def run_command(command):
    """Runs a command, logs its output in real-time, and returns success status."""
    command_str = ' '.join(map(str, command))
    print(f"{Colors.CYAN}Running Command:{Colors.ENDC} {command_str}")
    logging.getLogger().handlers[0].handle(logging.makeLogRecord({
        'msg': f"Running Command: {command_str}", 'levelname': 'INFO', 'levelno': logging.INFO,
        'pathname': '', 'lineno': 0, 'exc_info': None, 'func': None
    }))
    try:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1, universal_newlines=True
        )
        with process.stdout:
            for line in iter(process.stdout.readline, ''):
                print(f"{Colors.GREY}\t{line.strip()}{Colors.ENDC}")
                logging.getLogger().handlers[0].handle(logging.makeLogRecord({
                    'msg': line.strip(), 'levelname': 'INFO', 'levelno': logging.INFO,
                    'pathname': '', 'lineno': 0, 'exc_info': None, 'func': None
                }))
        return_code = process.wait()
        if return_code != 0:
            logging.error(f"Command failed with exit code {return_code}: {command_str}")
            return False
    except FileNotFoundError:
        logging.error(f"Error: The command '{command[0]}' was not found.")
        return False
    except Exception as e:
        logging.error(f"An unexpected error occurred while running command: {e}")
        return False
    logging.info("Command finished successfully.")
    return True

def file_exists_and_is_not_empty(file_path):
    """Checks if a file exists and is not empty."""
    return file_path.is_file() and file_path.stat().st_size > 0

def format_duration(seconds):
    """Formats a duration in seconds into a human-readable string."""
    if seconds < 60:
        return f"{seconds:.2f} seconds"
    
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    
    parts = []
    if hours > 0:
        parts.append(f"{int(hours)} hours")
    if minutes > 0:
        parts.append(f"{int(minutes)} minutes")
    if seconds > 0:
        parts.append(f"{int(seconds)} seconds")
        
    return ", ".join(parts)

def get_data_path(filename):
    """
    Finds the absolute path to a data file packaged with the installation.
    This is now compatible with Python 3.8, 3.9, and newer.
    """
    try:
        if sys.version_info < (3, 9):
            raise ImportError("Use pkg_resources for older Python versions")
        
        with resources.as_file(resources.files('amrscan').joinpath(filename)) as path:
            return path
    except (ImportError, AttributeError):
        return pkg_resources.resource_filename('amrscan', filename)

def main():
    start_time = time.monotonic()

    parser = argparse.ArgumentParser(
        description=f"amrscan pipeline v{__version__}. A wrapper script to run the amrscan analysis pipeline.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    # --- I/O Arguments ---
    parser.add_argument("-i", "--input-fastqs", required=True, type=str, help="Comma-separated list of input FASTQ files.")
    parser.add_argument("-o", "--output-prefix", required=True, help="Prefix for output files and the main output directory.")
    
    # --- Database Arguments (now optional) ---
    parser.add_argument("--db-card", default=None, help="Path to a custom CARD AMR database FASTA file (overrides packaged DB).")
    parser.add_argument("--db-scg", default=None, help="Path to a custom Single Copy Genes FASTA file (overrides packaged DB).")
    parser.add_argument("--db-scg-lengths", default=None, help="Path to a custom SCG gene lengths TSV file (overrides packaged DB).")
    parser.add_argument("--db-card-metadata", default=None, help="Path to a custom CARD metadata file (overrides packaged DB).")
    
    # --- Performance & Filtering Arguments ---
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads to use (default: 8).")
    parser.add_argument("--homscan-pid-cutoff", type=float, default=0.9, help="Minimum percent identity for HOMSCAN hits (0.0-1.0 scale). Default: 0.9")
    parser.add_argument("--varscan-pid-cutoff", type=float, default=0.9, help="Minimum nucleotide percent identity for VARSCAN hits (0.0-1.0 scale). Default: 0.9")
    parser.add_argument("--consensus-cutoff", type=float, default=0.9, help="Consensus cutoff for homscan family assignment. Default: 0.9")
    
    # --- Gene Type & PID Type Arguments ---
    parser.add_argument("--homscan-gene-types", default='H', help="Comma-delimited list of gene types for homscan (e.g., 'H,K'). Default: 'H'")
    parser.add_argument("--varscan-gene-types", default='V,R', help="Comma-delimited list of variant types for varscan (e.g., 'V,R,O'). Default: 'V,R'")
    parser.add_argument("--homscan-pid-type", choices=['protein', 'nucleotide'], default='protein', help="PID type to use for homscan filtering and WTA. Default: protein")
    
    # --- Control Flow Arguments ---
    parser.add_argument("--skip-mapping", action="store_true", help="Skip mapping and start from homscan/varscan.")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode, creating detailed logs for sub-scripts.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output directory.")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args()

    # --- Setup Paths ---
    fastq_files = [Path(f.strip()) for f in args.input_fastqs.split(',')]
    input_fastqs_str = ",".join(map(str, fastq_files))
    output_prefix_path = Path(args.output_prefix).resolve()
    output_prefix_name = output_prefix_path.name
    main_output_dir = output_prefix_path

    # --- Handle Existing Directory and --overwrite Flag ---
    if main_output_dir.exists():
        if not args.overwrite:
            print(f"{Colors.FAIL}Error: Output directory '{main_output_dir}' already exists.{Colors.ENDC}")
            print(f"{Colors.FAIL}Please remove it or use the --overwrite flag to proceed.{Colors.ENDC}")
            sys.exit(1)
        else:
            if not args.skip_mapping:
                print(f"{Colors.WARNING}Warning: --overwrite flag is set. The entire directory '{main_output_dir}' will be removed.{Colors.ENDC}")
                try:
                    shutil.rmtree(main_output_dir)
                except OSError as e:
                    print(f"{Colors.FAIL}Error: Could not remove directory {main_output_dir}. Reason: {e}{Colors.ENDC}")
                    sys.exit(1)
            else:
                print(f"{Colors.WARNING}Warning: --overwrite and --skip-mapping are set.{Colors.ENDC}")
                print(f"{Colors.WARNING}Cleaning up previous analysis results but preserving mapping data...{Colors.ENDC}")
                dirs_to_remove = [
                    main_output_dir / "tmp" / "homscan", main_output_dir / "tmp" / "varscan",
                    main_output_dir / "logs", main_output_dir / f"{output_prefix_name}_homscan_html",
                    main_output_dir / f"{output_prefix_name}_varscan_html"
                ]
                for d in dirs_to_remove:
                    if d.exists(): shutil.rmtree(d)
                for f in glob.glob(str(main_output_dir / f"{output_prefix_name}_*.tsv")):
                    os.remove(f)

    # --- Setup Logging and Directories ---
    log_dir = main_output_dir / "logs"
    setup_logging(log_dir / f"{output_prefix_name}_run.log")
    
    # --- Write Run Details ---
    run_details_file = log_dir / f"{output_prefix_name}_run_details.txt"
    print(f"{Colors.BLUE}{Colors.BOLD}--- amrscan pipeline starting (v{__version__}) ---{Colors.ENDC}")
    with open(run_details_file, "w") as f:
        f.write(f"AMRscan Pipeline Run Details\n" + "-"*30 + "\n")
        f.write(f"Pipeline Version: {__version__}\n")
        f.write(f"Date and Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Command-line: {' '.join(sys.argv)}\n\n")
        f.write("--- System & Tool Information ---\n")
        f.write(f"Python Version: {platform.python_version()}\n")
        for tool, version in get_tool_versions().items(): f.write(f"{tool} Version: {version}\n")
        for dep, version in get_py_dependency_versions().items(): f.write(f"{dep} Version: {version}\n")

    # --- Resolve Database Paths ---
    db_card = args.db_card or get_data_path('databases/amrscan_DB_CARD_v4.0.0/amrscan_DB_CARD_v4.0.0_NR_all_nuc.fasta')
    db_scg = args.db_scg or get_data_path('databases/amrscan_DB_SCG/SCG.fasta')
    db_scg_lengths = args.db_scg_lengths or get_data_path('databases/amrscan_DB_SCG/SCG_gene_lengths.tsv')
    db_card_metadata = args.db_card_metadata or get_data_path('databases/amrscan_DB_CARD_v4.0.0/amrscan_DB_CARD_v4.0.0_NR_metadata_curated.txt')

    # --- Define all sub-directories and script paths ---
    src_dir = Path(__file__).parent.resolve()
    tmp_dir = main_output_dir / "tmp"
    amrscan_tmp = tmp_dir / "amrscan"
    scgscan_tmp = tmp_dir / "scgscan"
    homscan_tmp = tmp_dir / "homscan"
    varscan_tmp = tmp_dir / "varscan"
    homscan_html_dir = main_output_dir / f"{output_prefix_name}_homscan_html"
    varscan_html_dir = main_output_dir / f"{output_prefix_name}_varscan_html"

    for d in [amrscan_tmp, scgscan_tmp, homscan_tmp, varscan_tmp, homscan_html_dir, varscan_html_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # ==========================================================================
    # MAPPING
    # ==========================================================================
    if not args.skip_mapping:
        log_boxed_header("Mapping", Colors.HEADER)
        amr_map_cmd = [sys.executable, str(src_dir / "amrscan_map_reads_bwa.py"), "-i", input_fastqs_str, "--db", db_card, "--tmp-dir", str(amrscan_tmp), "-t", str(args.threads)]
        if not run_command(amr_map_cmd): sys.exit(1)

        scg_map_cmd = [sys.executable, str(src_dir / "scgscan_map_reads_diamond.py"), "-i", input_fastqs_str, "--db", db_scg, "--tmp-dir", str(scgscan_tmp), "-t", str(args.threads)]
        if not run_command(scg_map_cmd): sys.exit(1)

        diamond_outputs = [scgscan_tmp / f"{f.name}.diamond.tsv" for f in fastq_files]
        uscg_report = scgscan_tmp / f"{output_prefix_name}_uscg_report.tsv"
        if all(file_exists_and_is_not_empty(p) for p in diamond_outputs):
            scg_quant_cmd = [sys.executable, str(src_dir / "scgscan_quantify_from_diamond.py"), "-i", ",".join(map(str, diamond_outputs)), "-l", db_scg_lengths, "-o", str(uscg_report), "-e", "1e-5"]
            if not run_command(scg_quant_cmd): sys.exit(1)
        else:
            logging.warning("One or more DIAMOND outputs for SCG not found or empty. Skipping SCG quantification.")
            uscg_report.touch()

        total_bases_file = scgscan_tmp / f"{output_prefix_name}_total_bases.txt"
        calc_bases_cmd = [sys.executable, str(src_dir / "scgscan_calculate_total_bases.py"), "-i", input_fastqs_str, "-o", str(total_bases_file)]
        if not run_command(calc_bases_cmd): sys.exit(1)
    else:
        print(f"{Colors.WARNING}Skipping mapping steps as requested by --skip-mapping.{Colors.ENDC}")
        total_bases_file = scgscan_tmp / f"{output_prefix_name}_total_bases.txt"
        uscg_report = scgscan_tmp / f"{output_prefix_name}_uscg_report.tsv"
        if not total_bases_file.is_file() or not uscg_report.is_file():
            logging.error("--skip-mapping was used, but required files from the mapping step are missing.")
            sys.exit(1)

    # --- Define dynamic file paths for downstream steps ---
    sam_files = [amrscan_tmp / f"{f.name}.sam" for f in fastq_files]
    sam_files_str = ",".join(map(str, sam_files))
    homscan_hits_files = [homscan_tmp / f"{s.name}_hits.tsv" for s in sam_files]
    homscan_hits_files_str = ",".join(map(str, homscan_hits_files))

    # ==========================================================================
    # HOMSCAN
    # ==========================================================================
    log_boxed_header("HomScan", Colors.HEADER)
    if all(file_exists_and_is_not_empty(p) for p in sam_files):
        homscan_process_cmd = [
            sys.executable, str(src_dir / "homscan_process_sam.py"), "-i", sam_files_str, "--db", db_card,
            "--tmp-dir", str(homscan_tmp), "--output-prefix", output_prefix_name,
            "--gene-types", args.homscan_gene_types
        ]
        if args.debug: homscan_process_cmd.append("--debug")
        if not run_command(homscan_process_cmd): sys.exit(1)

        if all(file_exists_and_is_not_empty(p) for p in homscan_hits_files):
            homscan_tabulate_cmd = [
                sys.executable, str(src_dir / "homscan_tabulate_and_normalise.py"), "-i", homscan_hits_files_str,
                "--metadata", db_card_metadata, "--uscg-report", str(uscg_report),
                "--total-bases-file", str(total_bases_file), "--pid-cutoff", str(args.homscan_pid_cutoff),
                "--consensus", str(args.consensus_cutoff), "--tmp-dir", str(homscan_tmp),
                "--output-prefix", output_prefix_name,
                "--pid-type", args.homscan_pid_type
            ]
            if not run_command(homscan_tabulate_cmd): sys.exit(1)

            homscan_resolve_cmd = [
                sys.executable, str(src_dir / "homscan_resolve_wta.py"), "-i", homscan_hits_files_str,
                "--metadata", db_card_metadata, "--pid-cutoff", str(args.homscan_pid_cutoff),
                "--tmp-dir", str(homscan_tmp), "--output-prefix", output_prefix_name,
                "--pid-type", args.homscan_pid_type
            ]
            if args.debug: homscan_resolve_cmd.append("--debug-wta")
            if not run_command(homscan_resolve_cmd): sys.exit(1)

            wta_summary = homscan_tmp / f"{output_prefix_name}_summary_wta.tsv"
            wta_assignments = homscan_tmp / f"{output_prefix_name}_wta_assignments.tsv"
            if file_exists_and_is_not_empty(wta_summary):
                homscan_visualise_cmd = [
                    sys.executable, str(src_dir / "homscan_visualise.py"), "-i", homscan_hits_files_str,
                    "--metadata", db_card_metadata, "--db", db_card,
                    "-o", str(homscan_html_dir), "--pid-cutoff", str(args.homscan_pid_cutoff),
                    "--wta-summary", str(wta_summary), "--wta-assignments", str(wta_assignments),
                    "--pid-type", args.homscan_pid_type
                ]
                if not run_command(homscan_visualise_cmd): sys.exit(1)
            else:
                logging.warning("WTA summary for homscan not found. Skipping visualization.")
        else:
            logging.warning("Homscan hits not found or empty. Skipping rest of homscan.")
    else:
        logging.warning("SAM files from AMR mapping not found or empty. Skipping homscan.")

    # ==========================================================================
    # VARSCAN
    # ==========================================================================
    log_boxed_header("VarScan", Colors.HEADER)
    if all(file_exists_and_is_not_empty(p) for p in sam_files):
        varscan_process_cmd = [
            sys.executable, str(src_dir / "varscan_process_sam.py"), "-i", sam_files_str, "--db", db_card,
            "--tmp-dir", str(varscan_tmp), "--output-prefix", output_prefix_name,
            "--gene-types", args.varscan_gene_types
        ]
        if args.debug: varscan_process_cmd.append("--debug")
        if not run_command(varscan_process_cmd): sys.exit(1)

        variant_hits = varscan_tmp / f"{output_prefix_name}_variant_hits.tsv"
        if file_exists_and_is_not_empty(variant_hits):
            varscan_tabulate_cmd = [
                sys.executable, str(src_dir / "varscan_tabulate_and_normalise.py"), "-i", str(variant_hits),
                "--metadata", db_card_metadata, "--uscg-report", str(uscg_report),
                "--total-bases-file", str(total_bases_file), "--tmp-dir", str(varscan_tmp),
                "--output-prefix", output_prefix_name, "--pid-cutoff", str(args.varscan_pid_cutoff)
            ]
            if not run_command(varscan_tabulate_cmd): sys.exit(1)

            variant_alignments = varscan_tmp / f"{output_prefix_name}_variant_alignments.txt"
            if file_exists_and_is_not_empty(variant_alignments):
                varscan_visualise_cmd = [
                    sys.executable, str(src_dir / "varscan_visualise.py"), "--variant-hits", str(variant_hits),
                    "--variant-alignments", str(variant_alignments), "--metadata", db_card_metadata,
                    "-o", str(varscan_html_dir)
                ]
                if not run_command(varscan_visualise_cmd): sys.exit(1)
            else:
                logging.warning("Variant alignments not found. Skipping varscan visualization.")
        else:
            logging.warning("Variant hits not found. Skipping rest of varscan.")
    else:
        logging.warning("SAM files from AMR mapping not found or empty. Skipping varscan.")

    # ==========================================================================
    # FINAL
    # ==========================================================================
    log_boxed_header("Final Report Consolidation", Colors.BLUE)
    final_cmd = [sys.executable, str(src_dir / "amrscan_consolidate_all.py"), "--output-prefix", str(output_prefix_path)]
    if not run_command(final_cmd): sys.exit(1)

    end_time = time.monotonic()
    duration_seconds = end_time - start_time
    
    print(f"\n{Colors.BLUE}{Colors.BOLD}--- amrscan pipeline finished successfully! ---{Colors.ENDC}")
    print(f"{Colors.BLUE}Total execution time: {format_duration(duration_seconds)}{Colors.ENDC}")

    print(f"\n{Colors.BOLD}--- Key Output Files ---{Colors.ENDC}")

    def summarize_files(title, directory, pattern):
        files = sorted(glob.glob(str(directory / pattern)))
        if files:
            print(f"\n{Colors.GREEN}{title}:{Colors.ENDC}")
            for f in files:
                print(f"  - {f}")
    
    summarize_files("Final Summary Tables", main_output_dir, f"{output_prefix_name}_*.tsv")
    summarize_files("HomScan Visualizations", homscan_html_dir, "*.html")
    summarize_files("VarScan Visualizations", varscan_html_dir, "*.html")
    summarize_files("Log Files", log_dir, "*")


if __name__ == "__main__":
    main()