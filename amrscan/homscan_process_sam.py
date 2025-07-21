#!/usr/bin/env python3

import argparse
import os
import sys
import re
import csv
from collections import defaultdict

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

# --- Start of Protein Translation Code ---

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '_', 'TAG': '_',
    'TGT': 'C', 'TGC': 'C', 'TGA': '_', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

def translate_codon(codon):
    """Translates a single codon. Handles gaps and ambiguous bases."""
    if len(codon) != 3 or '-' in codon or 'N' in codon.upper():
        return 'X'
    return CODON_TABLE.get(codon.upper(), 'X')

def calculate_protein_metrics(cigar, pos, query_seq, ref_nuc_seq):
    """
    Calculates protein PID, mismatches, and the denominator (codons compared).
    Returns a tuple: (pid, mismatches, denominator).
    """
    aligned_ref_nuc, aligned_query_nuc = [], []
    query_cursor, ref_pos = 0, pos - 1

    cigar_parts = re.findall(r'(\d+)([MDN=XISH])', cigar)
    for length, op in cigar_parts:
        length = int(length)
        if op in ('M', '=', 'X'):
            aligned_ref_nuc.append(ref_nuc_seq[ref_pos : ref_pos + length])
            aligned_query_nuc.append(query_seq[query_cursor : query_cursor + length])
            ref_pos += length
            query_cursor += length
        elif op == 'I':
            aligned_ref_nuc.append('-' * length)
            aligned_query_nuc.append(query_seq[query_cursor : query_cursor + length])
            query_cursor += length
        elif op == 'D':
            aligned_ref_nuc.append(ref_nuc_seq[ref_pos : ref_pos + length])
            aligned_query_nuc.append('-' * length)
            ref_pos += length
        elif op in ('S', 'H'):
            query_cursor += length

    frame_offset = (pos - 1) % 3
    padded_ref_str = ('-' * frame_offset) + "".join(aligned_ref_nuc)
    padded_query_str = ('-' * frame_offset) + "".join(aligned_query_nuc)

    if len(padded_ref_str) % 3 != 0:
        padding_needed = 3 - (len(padded_ref_str) % 3)
        padded_ref_str += '-' * padding_needed
        padded_query_str += '-' * padding_needed

    matches, codons_compared = 0, 0
    for i in range(0, len(padded_ref_str), 3):
        ref_codon = padded_ref_str[i:i+3]
        
        if len(ref_codon) == 3 and '-' not in ref_codon:
            ref_aa = translate_codon(ref_codon)
            
            if ref_aa not in ('_', 'X'):
                codons_compared += 1
                query_codon = padded_query_str[i:i+3]
                query_aa = translate_codon(query_codon)
                if ref_aa == query_aa:
                    matches += 1
    
    if codons_compared == 0:
        return (0.0, 0, 0)
    
    pid = (matches / codons_compared) * 100.0
    mismatches = codons_compared - matches
    return (pid, mismatches, codons_compared)

def parse_fasta_db(db_path):
    """Parses a FASTA file to get sequence lengths and nucleotide sequences."""
    print(BColors.cyan(f"--- Parsing reference database from: {db_path} ---"))
    if not os.path.exists(db_path):
        print(BColors.red(f"Error: Database file not found at '{db_path}'"), file=sys.stderr)
        sys.exit(1)
        
    lengths, sequences = {}, {}
    current_id, current_seq_parts = None, []

    with open(db_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    nuc_seq = "".join(current_seq_parts)
                    sequences[current_id] = nuc_seq
                    lengths[current_id] = len(nuc_seq)
                current_id = line[1:] 
                current_seq_parts = []
            else:
                current_seq_parts.append(line.upper())
    
    if current_id:
        nuc_seq = "".join(current_seq_parts)
        sequences[current_id] = nuc_seq
        lengths[current_id] = len(nuc_seq)
        
    print(BColors.green(f"--- Found and processed {len(sequences)} sequences in the database."))
    return lengths, sequences

def get_nm_tag(optional_fields):
    for field in optional_fields:
        if field.startswith('NM:i:'):
            try: return int(field[5:])
            except (ValueError, IndexError): return None
    return None
    
def get_alignment_length_from_cigar(cigar):
    parts = re.findall(r'(\d+)([MDN=X])', cigar)
    return sum(int(num) for num, op in parts)

def is_valid_alignment(cigar, pos, aln_len, target_len):
    """
    Checks if an alignment is valid under the strict full-read alignment policy.
    Clipping is only allowed if the alignment reaches the end of the reference.
    """
    has_clipping = 'S' in cigar or 'H' in cigar
    if not has_clipping:
        return True

    cigar_ops = re.findall(r'(\d+)([A-Z])', cigar)
    
    if cigar_ops[0][1] in ('S', 'H') and len(cigar_ops) == 2:
        if pos == 1:
            return True

    if cigar_ops[-1][1] in ('S', 'H') and len(cigar_ops) == 2:
        if (pos + aln_len - 1) == target_len:
            return True
            
    if cigar_ops[0][1] in ('S', 'H') and cigar_ops[-1][1] in ('S', 'H'):
        if pos == 1 and (pos + aln_len - 1) == target_len:
            return True

    return False

def process_sam_file(sam_path, db_lengths, db_sequences, tmp_dir, min_aln_len, file_idx, debug_writer, allowed_gene_types):
    """
    Processes a single SAM file for AMR hits, appending a file-specific suffix to query_ids.
    """
    print(BColors.cyan(f"\n--- Processing SAM file for AMR hits: {os.path.basename(sam_path)} (Input File #{file_idx}) ---"))
    
    query_sequences = {}

    with open(sam_path, 'r') as f:
        for line in f:
            if line.startswith('@'): continue
            try:
                fields = line.strip().split('\t')
                original_qname = fields[0]
                suffixed_qname = f"{original_qname}_{file_idx}"
                seq = fields[9]
                if seq != '*':
                    query_sequences[suffixed_qname] = seq
            except IndexError: continue

    input_basename = os.path.basename(sam_path)
    hits_path = os.path.join(tmp_dir, f"{input_basename}_hits.tsv")
    
    hits_by_query_id = defaultdict(int)
    with open(hits_path, 'w') as f_hits:
        header = [
            "ARO_matched", "query_id",
            "nucleotide_pid", "nucleotide_mismatches", "nucleotide_denominator",
            "protein_pid", "protein_mismatches", "protein_denominator",
            "position_on_ref", "read_length", "target_length"
        ]
        f_hits.write('\t'.join(header) + '\n')

        with open(sam_path, 'r') as f_in:
            for line in f_in:
                if line.startswith('@'): continue
                
                try:
                    fields = line.strip().split('\t')
                    original_qname = fields[0]
                    suffixed_qname = f"{original_qname}_{file_idx}"
                    
                    actual_seq = query_sequences.get(suffixed_qname)
                    if not actual_seq: continue

                    rname = fields[2]
                    if rname == '*': continue

                    gene_type_match = re.search(r'__([HKVORU])(?:__|$)', rname)
                    if not gene_type_match:
                        print(BColors.yellow(f"Warning: Skipping SAM line due to unparseable gene type in RNAME: {rname}"), file=sys.stderr)
                        continue
                    gene_type = gene_type_match.group(1)

                    debug_row = None
                    if debug_writer:
                        debug_row = {
                            'query_id': suffixed_qname, 'ARO_matched': rname,
                            'status': 'N/A', 'reason': 'N/A', 'details': 'N/A',
                            'nucleotide_pid': 'N/A', 'aln_len': 'N/A'
                        }

                    if gene_type in allowed_gene_types:
                        pos_str = fields[3]
                        cigar = fields[5]
                        optional_fields = fields[11:]
                        
                        aln_len = get_alignment_length_from_cigar(cigar)
                        if debug_writer:
                            debug_row['aln_len'] = aln_len

                        if aln_len < min_aln_len:
                            if debug_writer:
                                debug_row['status'] = 'Rejected'
                                debug_row['reason'] = 'Alignment length below cutoff'
                                debug_row['details'] = f"Aln length {aln_len} < {min_aln_len}"
                                debug_writer.writerow(debug_row)
                            continue

                        pos = int(pos_str)
                        target_length = db_lengths.get(rname)
                        
                        if not is_valid_alignment(cigar, pos, aln_len, target_length):
                            if debug_writer:
                                debug_row['status'] = 'Rejected'
                                debug_row['reason'] = 'Does not meet full read alignment criteria'
                                debug_row['details'] = f"CIGAR: {cigar}, POS: {pos}, REF_LEN: {target_length}"
                                debug_writer.writerow(debug_row)
                            continue

                        read_length = len(actual_seq)
                        
                        num_mismatches = get_nm_tag(optional_fields)
                        nuc_pid_str = "NA"
                        if num_mismatches is not None and aln_len > 0:
                            pid_val = ((aln_len - num_mismatches) / aln_len) * 100
                            nuc_pid_str = f"{pid_val:.2f}"
                            if debug_writer:
                                debug_row['nucleotide_pid'] = nuc_pid_str
                        
                        nuc_mismatches_str = str(num_mismatches) if num_mismatches is not None else "NA"

                        prot_pid_str, prot_mismatches_str, prot_denom_str = "NA", "NA", "NA"
                        ref_nuc_seq = db_sequences.get(rname)
                        if ref_nuc_seq:
                            prot_pid, prot_mismatches, prot_denom = calculate_protein_metrics(cigar, pos, actual_seq, ref_nuc_seq)
                            prot_pid_str = f"{prot_pid:.2f}"
                            prot_mismatches_str = str(prot_mismatches)
                            prot_denom_str = str(prot_denom)
                            
                            try:
                                nuc_pid = float(nuc_pid_str)
                                prot_pid = float(prot_pid_str)
                                if prot_pid < nuc_pid:
                                    prot_pid_str = nuc_pid_str
                            except ValueError:
                                pass
                        
                        if debug_writer:
                            debug_row['status'] = 'Accepted'
                            debug_row['reason'] = 'Passes all filters'
                            debug_writer.writerow(debug_row)

                        output_data = [
                            rname, suffixed_qname,
                            nuc_pid_str, nuc_mismatches_str, str(aln_len),
                            prot_pid_str, prot_mismatches_str, prot_denom_str,
                            pos_str, read_length, str(target_length or "NA")
                        ]
                        f_hits.write('\t'.join(map(str, output_data)) + '\n')
                        hits_by_query_id[suffixed_qname] += 1
                    
                    elif debug_writer:
                        debug_row['status'] = 'Skipped'
                        debug_row['reason'] = 'Gene type not in specified list'
                        debug_row['details'] = f"Gene type was '{gene_type}', allowed types are {list(allowed_gene_types)}"
                        debug_writer.writerow(debug_row)


                except (IndexError, ValueError) as e:
                    print(BColors.yellow(f"Warning: Skipping malformed SAM line: {line.strip()} | Error: {e}"), file=sys.stderr)
                    continue
    
    print(BColors.cyan("--- Adding uniqueness information to the hits table ---"))
    try:
        with open(hits_path, 'r') as f_in:
            lines = f_in.readlines()
        
        with open(hits_path, 'w') as f_out:
            header = lines[0].strip()
            f_out.write(header + "\tuniqueness\n")
            for line in lines[1:]:
                fields = line.strip().split('\t')
                query_id = fields[1]
                uniqueness = "UNIQUE" if hits_by_query_id[query_id] == 1 else "NOT_UNIQUE"
                f_out.write('\t'.join(fields) + f"\t{uniqueness}\n")
    except Exception as e:
        print(BColors.red(f"Error adding uniqueness information to {hits_path}: {e}"), file=sys.stderr)
        sys.exit(1)
    
    files_created = [hits_path]
    print(BColors.green("--- Analysis complete. Generated files: ---"))
    for f in files_created: print(f"- {f}")
    return files_created

def main():
    parser = argparse.ArgumentParser(description="Analyzes SAM files for AMR hits against an AMR database.")
    parser.add_argument("-i", "--input-files", required=True, help="Comma-delimited list of input SAM files.")
    parser.add_argument("-d", "--db", required=True, help="Path to the reference database FASTA file.")
    parser.add_argument("--tmp-dir", default=".", help="Directory for output files.")
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix for the output debug file (e.g., 'MyProject')."
    )
    parser.add_argument(
        "--min-aln-len",
        type=int,
        default=100,
        help="Minimum alignment length (in base pairs) for a hit to be reported. Default: 100"
    )
    parser.add_argument(
        "--gene-types",
        default='H',
        help="Comma-delimited list of gene types to process (e.g., 'H,K'). Default: 'H'"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug mode, which creates a detailed log file of all homology checks."
    )
    args = parser.parse_args()

    os.makedirs(args.tmp_dir, exist_ok=True)
    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    if not input_files_list:
        print(BColors.red("Error: No input SAM files provided."), file=sys.stderr)
        sys.exit(1)

    db_lengths, db_sequences = parse_fasta_db(args.db)

    allowed_gene_types = {t.strip().upper() for t in args.gene_types.split(',') if t.strip()}
    print(BColors.cyan(f"--- Processing gene types: {', '.join(sorted(list(allowed_gene_types)))} ---"))

    debug_writer = None
    debug_file = None
    if args.debug:
        debug_output_path = os.path.join(args.tmp_dir, f"{args.output_prefix}_homology_debug.tsv")
        print(BColors.yellow(f"--- DEBUG MODE ON. Writing detailed log to: {debug_output_path} ---"))
        debug_file = open(debug_output_path, 'w', newline='')
        header = ['query_id', 'ARO_matched', 'status', 'reason', 'details', 'nucleotide_pid', 'aln_len']
        debug_writer = csv.DictWriter(debug_file, fieldnames=header, delimiter='\t')
        debug_writer.writeheader()

    all_created_files = []
    for i, sam_file in enumerate(input_files_list):
        if not os.path.exists(sam_file):
            print(BColors.yellow(f"Warning: Input file not found, skipping: {sam_file}"), file=sys.stderr)
            continue
        created = process_sam_file(sam_file, db_lengths, db_sequences, args.tmp_dir, args.min_aln_len, i + 1, debug_writer, allowed_gene_types)
        all_created_files.extend(created)

    if debug_file:
        debug_file.close()

    print(BColors.green("\n\n--- All Processing Complete ---"))
    print("The following output files were generated in total:")
    for f in all_created_files: print(f"- {f}")

if __name__ == "__main__":
    main()