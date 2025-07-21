#!/usr/bin/env python3

import requests
import gzip
import os
import shutil
import argparse

# <<< THIS IS THE ONLY FUNCTION THAT HAS BEEN MODIFIED >>>
def process_fasta_entry(source_fasta_path, new_identifier, target_fasta_path):
    """
    Reads a source FASTA file, modifies sequence headers by appending a unique,
    padded index to the new_identifier, and saves to a target FASTA file.

    Args:
        source_fasta_path (str): Path to the input FASTA file (e.g., COG0001.fa).
        new_identifier (str): The new identifier to use as a base (e.g., argS).
        target_fasta_path (str): Path to save the modified FASTA file (e.g., argS.fa).

    Returns:
        bool: True if processing was successful, False otherwise.
    """
    try:
        with open(source_fasta_path, 'r') as infile, open(target_fasta_path, 'w') as outfile:
            print(f"Processing '{os.path.basename(source_fasta_path)}' with ID '{new_identifier}' -> '{os.path.basename(target_fasta_path)}'")
            
            index = 1 # Initialize a counter for each file
            for line in infile:
                if line.startswith(">"):
                    # Format the index with leading zeros to 3 digits (001, 002, etc.)
                    padded_index = f"{index:03d}"
                    original_header_content = line[1:].strip()
                    
                    # Construct the new header with the padded index
                    # Note the change from a tab to a space after the new ID
                    outfile.write(f">{new_identifier}___{padded_index} {original_header_content}\n")
                    index += 1 # Increment for the next sequence
                else:
                    outfile.write(line)
        print(f"Successfully created '{os.path.basename(target_fasta_path)}'")
        return True
    except IOError as e:
        print(f"Error processing file '{os.path.basename(source_fasta_path)}' or writing to '{os.path.basename(target_fasta_path)}': {e}")
        return False
    except Exception as e:
        print(f"An unexpected error occurred during FASTA processing for '{os.path.basename(source_fasta_path)}': {e}")
        if os.path.exists(target_fasta_path):
            try:
                os.remove(target_fasta_path)
            except OSError:
                pass
        return False

def download_and_process_cog_fastas(cog_name_pairs, output_dir):
    """
    For each (COG_ID, NewName) pair:
    1. Downloads the FASTA file for COG_ID from NCBI (if not already present).
    2. Decompresses it to {COG_ID}.fa, deleting the .gz file.
    3. Modifies sequence headers in {COG_ID}.fa using NewName.
    4. Saves the result as {NewName}.fa in the output_dir.
    Intermediate {COG_ID}.fa files are kept.
    5. After all processing, intermediate {COG_ID}.fa files are deleted unless they
       were also specified as a final target name.
    Args:
        cog_name_pairs (list): A list of tuples (cog_id, new_name).
        output_dir (str): The directory for all FASTA files (intermediate and final).
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {os.path.abspath(output_dir)}")
    else:
        print(f"Output directory: {os.path.abspath(output_dir)}")

    base_url = "https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fasta/{}.fa.gz"
    success_count = 0
    error_count = 0
    downloaded_intermediate_cogs = set()

    for cog_id, new_name in cog_name_pairs:
        if not cog_id or not new_name:
            print(f"\nSkipping invalid entry: COG ID '{cog_id}', New Name '{new_name}'")
            error_count +=1
            continue

        file_url = base_url.format(cog_id)

        gzipped_cog_filename = f"{cog_id}.fa.gz"
        intermediate_cog_filename = f"{cog_id}.fa"
        final_target_filename = f"{new_name}.fa"

        gz_filepath = os.path.join(output_dir, gzipped_cog_filename)
        intermediate_cog_filepath = os.path.join(output_dir, intermediate_cog_filename)
        final_target_filepath = os.path.join(output_dir, final_target_filename)

        print(f"\nProcessing COG ID '{cog_id}' for target '{final_target_filename}'...")

        if os.path.exists(final_target_filepath):
            print(f"Final target '{final_target_filename}' already exists in '{os.path.abspath(output_dir)}'. Skipping.")
            success_count += 1
            continue

        intermediate_ready = False
        if cog_id in downloaded_intermediate_cogs:
            if os.path.exists(intermediate_cog_filepath):
                print(f"Using previously downloaded/decompressed intermediate '{intermediate_cog_filename}'.")
                intermediate_ready = True
            else:
                print(f"Error: '{cog_id}' was marked as downloaded, but '{intermediate_cog_filename}' not found. Attempting re-download.")
                downloaded_intermediate_cogs.discard(cog_id)
        elif os.path.exists(intermediate_cog_filepath):
             print(f"Found existing intermediate '{intermediate_cog_filename}'.")
             intermediate_ready = True
             if os.path.exists(gz_filepath):
                try:
                    os.remove(gz_filepath)
                    print(f"Removed lingering '{gzipped_cog_filename}'.")
                except OSError as e:
                    print(f"Could not remove lingering '{gzipped_cog_filename}': {e}")

        if not intermediate_ready:
            print(f"Intermediate '{intermediate_cog_filename}' not found. Attempting download...")
            try:
                print(f"Downloading '{gzipped_cog_filename}' from {file_url}...")
                response = requests.get(file_url, stream=True, timeout=60)
                response.raise_for_status()

                with open(gz_filepath, 'wb') as f_gz:
                    shutil.copyfileobj(response.raw, f_gz)
                print(f"Successfully downloaded '{gzipped_cog_filename}'.")

                print(f"Decompressing '{gzipped_cog_filename}' to '{intermediate_cog_filename}'...")
                with gzip.open(gz_filepath, 'rb') as f_in:
                    with open(intermediate_cog_filepath, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                print(f"Successfully decompressed to '{intermediate_cog_filename}'.")
                downloaded_intermediate_cogs.add(cog_id)
                intermediate_ready = True

            except requests.exceptions.HTTPError as e:
                print(f"HTTP error downloading {cog_id}: {e}")
                if e.response.status_code == 404:
                    print(f"This may mean COG ID '{cog_id}' does not exist at the source URL or the URL is incorrect.")
                error_count += 1
            except requests.exceptions.RequestException as e:
                print(f"Error downloading {cog_id}: {e}")
                error_count += 1
            except gzip.BadGzipFile:
                print(f"Error: '{gzipped_cog_filename}' is not a valid GZIP file or is corrupted.")
                print(f"Please try deleting '{gz_filepath}' if it exists and running again.")
                error_count += 1
            except Exception as e:
                print(f"An unexpected error occurred during download/decompression of {cog_id}: {e}")
                error_count += 1
            finally:
                if os.path.exists(gz_filepath):
                    try:
                        os.remove(gz_filepath)
                        if intermediate_ready:
                            print(f"Removed '{gzipped_cog_filename}'.")
                        else:
                            print(f"Cleaned up '{gzipped_cog_filename}' after error.")
                    except OSError as e_rm:
                        print(f"Error removing '{gzipped_cog_filename}': {e_rm}")

                if not intermediate_ready and os.path.exists(intermediate_cog_filepath):
                    try:
                        os.remove(intermediate_cog_filepath)
                        print(f"Cleaned up partially created '{intermediate_cog_filename}' after error.")
                    except OSError:
                        pass

        if intermediate_ready:
            if process_fasta_entry(intermediate_cog_filepath, new_name, final_target_filepath):
                success_count += 1
            else:
                error_count += 1
        elif not os.path.exists(final_target_filepath):
            print(f"Skipping processing for '{final_target_filename}' due to earlier errors with '{cog_id}'.")

    print(f"\n--- Cleaning up intermediate COG FASTA files ---")
    final_target_filenames_set = {f"{pair[1]}.fa" for pair in cog_name_pairs}
    all_input_cog_ids = {pair[0] for pair in cog_name_pairs}

    cleaned_intermediate_count = 0
    cleanup_error_count = 0

    for cog_id_to_clean in all_input_cog_ids:
        intermediate_filename = f"{cog_id_to_clean}.fa"
        intermediate_filepath = os.path.join(output_dir, intermediate_filename)

        if intermediate_filename not in final_target_filenames_set:
            if os.path.exists(intermediate_filepath):
                try:
                    os.remove(intermediate_filepath)
                    print(f"Successfully removed intermediate file: {intermediate_filepath}")
                    cleaned_intermediate_count += 1
                except OSError as e:
                    print(f"Error removing intermediate file {intermediate_filepath}: {e}")
                    cleanup_error_count += 1

    if len(all_input_cog_ids) > 0:
        print(f"Intermediate file cleanup summary: {cleaned_intermediate_count} removed, {cleanup_error_count} errors.")

    print(f"\n--- Summary ---")
    print(f"Total COG ID/New Name pairs requested: {len(cog_name_pairs)}")
    print(f"Successfully created/found existing final FASTA files: {success_count}")
    print(f"Errors encountered for COG ID/New Name pairs: {error_count}")
    if error_count > 0:
        print("Please check the error messages above for details on failed downloads/processing.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="amrscan_DB_USCG_download.py",
        description=(
            "Downloads COG FASTA files from NCBI, decompresses them, "
            "modifies sequence headers based on a two-column list file, "
            "and saves them with new names. Intermediate {COG_ID}.fa files are then deleted."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-d", "--directory",
        default="cog_fastas",
        help=("The directory to save the final processed FASTA files ({NewName}.fa).\n"
              "Intermediate {COG_ID}.fa files are temporarily stored here during processing and then deleted.\n"
              "Defaults to 'cog_fastas'.")
    )
    parser.add_argument(
        "-l", "--list-file",
        required=True,
        help=(
            "Path to a tab-separated file containing COG IDs and their desired output names.\n"
            "Each line should have two columns: COG_ID\\tNewName (e.g., 'COG0001\\tMyGene1').\n"
            "This argument is required."
        )
    )
    args = parser.parse_args()

    cog_name_pairs_to_process = []

    try:
        with open(args.list_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split('\t')
                if len(parts) == 2:
                    cog_id = parts[0].strip()
                    new_name = parts[1].strip()
                    if cog_id and new_name:
                        cog_name_pairs_to_process.append((cog_id, new_name))
                    else:
                        print(f"Warning: Skipping line {line_num} in '{args.list_file}' due to empty COG ID or New Name after stripping: '{line}'")
                else:
                    print(f"Warning: Skipping malformed line {line_num} in '{args.list_file}' (expected 2 tab-separated columns): '{line}'")

        if not cog_name_pairs_to_process:
            print(f"Error: The list file '{args.list_file}' is empty or contains no valid COG ID/New Name pairs.")
            exit(1)
        print(f"Loaded {len(cog_name_pairs_to_process)} COG ID/New Name pair(s) from '{args.list_file}'.")

    except FileNotFoundError:
        print(f"Error: The specified list file '{args.list_file}' was not found.")
        exit(1)
    except Exception as e:
        print(f"Error reading COG list file '{args.list_file}': {e}")
        exit(1)

    download_and_process_cog_fastas(cog_name_pairs_to_process, args.directory)