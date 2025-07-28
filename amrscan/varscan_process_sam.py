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

# --- Protein Translation Code ---
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

def translate_sequence(nuc_sequence):
    """Translates a nucleotide sequence into an amino acid sequence."""
    nuc_sequence = nuc_sequence.replace('-', '')
    aa_sequence = []
    for i in range(0, len(nuc_sequence) - (len(nuc_sequence) % 3), 3):
        aa_sequence.append(translate_codon(nuc_sequence[i:i+3]))
    return "".join(aa_sequence)

def parse_fasta_db(db_path):
    """Parses a FASTA file to get reference sequences."""
    print(BColors.cyan(f"--- Parsing reference database from: {db_path} ---"))
    if not os.path.exists(db_path):
        print(BColors.red(f"Error: Database file not found at '{db_path}'"), file=sys.stderr)
        sys.exit(1)
        
    sequences = {}
    current_id, current_seq_parts = None, []

    with open(db_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = "".join(current_seq_parts)
                current_id = line[1:]
                current_seq_parts = []
            else:
                current_seq_parts.append(line.upper())
    
    if current_id:
        sequences[current_id] = "".join(current_seq_parts)
        
    print(BColors.green(f"--- Found and processed {len(sequences)} sequences in the database."))
    return sequences

def parse_mutation_info_from_db(db_sequences):
    """Parses mutation details from reference sequence headers."""
    print(BColors.cyan("--- Parsing mutation information from database headers ---"))
    mutation_db = {}
    mutation_pattern = re.compile(r"([A-Z])(\d+)([A-Z*])", re.IGNORECASE)

    for header in db_sequences.keys():
        parts = header.split('__')
        if len(parts) < 4: continue
        
        aro_id = header
        gene_type = parts[-2]
        mutation_str = parts[-1]

        if gene_type in {'O', 'V', 'R'}:
            match = mutation_pattern.match(mutation_str)
            if match:
                ref_char, pos_str, mut_char = match.groups()
                mutation_db[aro_id] = {
                    'type': gene_type,
                    'pos': int(pos_str),
                    'ref': ref_char,
                    'mut': mut_char,
                    'full_mutation': mutation_str
                }
    print(BColors.green(f"--- Found mutation information for {len(mutation_db)} variant AROs."))
    return mutation_db

def get_nm_tag(optional_fields):
    """Extracts the NM (edit distance) tag from SAM optional fields."""
    for field in optional_fields:
        if field.startswith('NM:i:'):
            try: return int(field[5:])
            except (ValueError, IndexError): return None
    return None

def get_reference_span_from_cigar(cigar):
    """Calculates alignment span on the reference from CIGAR string (M, D, N, X, =)."""
    parts = re.findall(r'(\d+)([MDN=X])', cigar)
    return sum(int(num) for num, op in parts)

def is_valid_alignment(cigar, pos, aln_len, target_len):
    """
    Checks if an alignment is valid under the strict full-read alignment policy.
    """
    has_clipping = 'S' in cigar or 'H' in cigar
    if not has_clipping:
        return True

    cigar_ops = re.findall(r'(\d+)([A-Z])', cigar)
    
    if cigar_ops[0][1] in ('S', 'H') and len(cigar_ops) == 2:
        if (pos + aln_len - 1) == target_len:
            return True

    if cigar_ops[-1][1] in ('S', 'H') and len(cigar_ops) == 2:
        if pos == 1:
            return True
            
    if cigar_ops[0][1] in ('S', 'H') and cigar_ops[-1][1] in ('S', 'H'):
        if pos == 1 and (pos + aln_len - 1) == target_len:
            return True

    return False

def check_mutation_presence(cigar, pos, query_seq, ref_seq, mutation_info):
    """
    Checks for a specific mutation and returns a detailed result dictionary.
    """
    result = { 'confirmed': False, 'reason': 'Unknown Error', 'details': 'N/A' }

    aligned_query, aligned_ref = "", ""
    query_cursor, ref_cursor = 0, pos - 1
    cigar_parts = re.findall(r'(\d+)([MDN=XISH])', cigar)

    for length, op in cigar_parts:
        length = int(length)
        if op in ('M', '=', 'X'):
            aligned_query += query_seq[query_cursor : query_cursor + length]
            aligned_ref += ref_seq[ref_cursor : ref_cursor + length]
            query_cursor += length
            ref_cursor += length
        elif op == 'I':
            aligned_query += query_seq[query_cursor : query_cursor + length]
            aligned_ref += '-' * length
            query_cursor += length
        elif op == 'D':
            aligned_query += '-' * length
            aligned_ref += ref_seq[ref_cursor : ref_cursor + length]
            ref_cursor += length
        elif op in ('S', 'H'):
            query_cursor += length

    mut_type = mutation_info['type']
    mut_pos_1based = mutation_info['pos']
    expected_mut_char = mutation_info['mut']

    if mut_type == 'R':
        target_ref_nuc_pos_0based = mut_pos_1based - 1
        ref_nuc_count = 0
        alignment_idx = -1
        for i, char in enumerate(aligned_ref):
            if char != '-':
                if (pos - 1 + ref_nuc_count) == target_ref_nuc_pos_0based:
                    alignment_idx = i
                    break
                ref_nuc_count += 1
        
        if alignment_idx == -1:
            result['reason'] = 'Rejected: Read does not cover nucleotide mutation site'
            return result
        
        query_nuc = aligned_query[alignment_idx]
        if query_nuc.upper() == expected_mut_char.upper():
            result['confirmed'] = True
            result['reason'] = 'Accepted: Confirmed nucleotide mutation'
            result['details'] = f"Found '{query_nuc}' at position {mut_pos_1based}"
        else:
            result['reason'] = 'Rejected: Nucleotide mismatch at mutation site'
            result['details'] = f"Found '{query_nuc}', expected '{expected_mut_char}' at position {mut_pos_1based}"
        return result

    elif mut_type in {'V', 'O'}:
        target_codon_start_0based = (mut_pos_1based - 1) * 3
        ref_nuc_count = 0
        codon_start_idx = -1
        for i, char in enumerate(aligned_ref):
            if char != '-':
                if (pos - 1 + ref_nuc_count) == target_codon_start_0based:
                    codon_start_idx = i
                    break
                ref_nuc_count += 1
        
        if codon_start_idx == -1 or (codon_start_idx + 3) > len(aligned_query):
            result['reason'] = 'Rejected: Read does not cover full codon site'
            return result
            
        query_codon_str = aligned_query[codon_start_idx : codon_start_idx + 3].replace('-', '')
        if len(query_codon_str) != 3:
            result['reason'] = 'Rejected: Indel within the target codon'
            result['details'] = f"Query codon '{query_codon_str}' is not 3bp long"
            return result

        query_aa = translate_codon(query_codon_str)
        if query_aa.upper() == expected_mut_char.upper():
            result['confirmed'] = True
            result['reason'] = 'Accepted: Confirmed amino acid mutation'
            result['details'] = f"Found '{query_aa}' at AA position {mut_pos_1based}"
        else:
            result['reason'] = 'Rejected: Amino acid mismatch at mutation site'
            result['details'] = f"Found '{query_aa}', expected '{expected_mut_char}' at AA position {mut_pos_1based}"
        return result

    result['reason'] = f"Rejected: Unknown mutation type '{mut_type}'"
    return result

def format_fullseq_alignment_for_output(cigar, full_query_seq, full_ref_seq, mutation_info, sam_pos, query_id, aro_id):
    """
    Creates a human-readable, gapped alignment showing the full reference and query 
    sequences with a marker pointing to the mutation site.
    """
    mut_type = mutation_info['type']
    mut_pos_1based = mutation_info['pos']
    if mut_type == 'R':
        target_ref_pos_0based = mut_pos_1based - 1
    else:
        target_ref_pos_0based = (mut_pos_1based - 1) * 3

    gapped_ref_block = ""
    gapped_que_block = ""
    marker_col_index = -1

    query_cursor = 0
    ref_cursor = sam_pos - 1
    alignment_col = 0

    cigar_parts = re.findall(r'(\d+)([MDN=XISH])', cigar)

    for length, op in cigar_parts:
        length = int(length)
        if op in ('M', '=', 'X'):
            for _ in range(length):
                if ref_cursor == target_ref_pos_0based and marker_col_index == -1:
                    marker_col_index = alignment_col
                gapped_ref_block += full_ref_seq[ref_cursor]
                gapped_que_block += full_query_seq[query_cursor]
                ref_cursor += 1
                query_cursor += 1
                alignment_col += 1
        elif op == 'I':
            for _ in range(length):
                gapped_ref_block += '-'
                gapped_que_block += full_query_seq[query_cursor]
                query_cursor += 1
                alignment_col += 1
        elif op == 'D':
            for _ in range(length):
                if ref_cursor == target_ref_pos_0based and marker_col_index == -1:
                    marker_col_index = alignment_col
                gapped_ref_block += full_ref_seq[ref_cursor]
                gapped_que_block += '-'
                ref_cursor += 1
                alignment_col += 1
        elif op == 'S':
            query_cursor += length

    ref_prefix = full_ref_seq[:sam_pos - 1]
    que_prefix = '.' * len(ref_prefix)
    
    ref_suffix = full_ref_seq[ref_cursor:]
    que_suffix = '.' * len(ref_suffix)

    final_ref_str = ref_prefix + gapped_ref_block + ref_suffix
    final_que_str = que_prefix + gapped_que_block + que_suffix

    output_lines = []
    output_lines.append(f"> {query_id} | maps to | {aro_id} | mutation: {mutation_info['full_mutation']}")
    output_lines.append(f"REFSEQ:    {final_ref_str}")
    output_lines.append(f"QUESEQ:    {final_que_str}")
    
    if marker_col_index != -1:
        final_marker_pos = len(ref_prefix) + marker_col_index
        marker_line = "MARKER:    " + " " * final_marker_pos + "^"
        output_lines.append(marker_line)
        
    output_lines.append("")
    
    return "\n".join(output_lines)

def process_sam_files_for_variants(sam_files, db_sequences, mutation_db, output_path, debug_mode, debug_output_path, alignment_output_path, min_aln_len, require_full_read_aln, allowed_variant_types):
    """
    Processes SAM files, filters by alignment quality, checks for specific mutations,
    and writes confirmed hits and debug files.
    """
    print(BColors.cyan(f"\n--- Processing {len(sam_files)} SAM file(s) for confirmed variants ---"))
    
    confirmed_hits = []
    
    debug_writer = None
    debug_file = None
    
    print(BColors.cyan(f"--- Writing confirmed hit alignments to: {alignment_output_path} ---"))
    align_file = open(alignment_output_path, 'w')

    if debug_mode:
        print(BColors.yellow(f"--- DEBUG MODE ON. Writing detailed log to: {debug_output_path} ---"))
        debug_file = open(debug_output_path, 'w', newline='')
        debug_header = [
            'query_id', 'ARO_matched', 'status', 'reason', 'details', 'nucleotide_pid'
        ]
        debug_writer = csv.DictWriter(debug_file, fieldnames=debug_header, delimiter='\t')
        debug_writer.writeheader()

    for file_idx, sam_path in enumerate(sam_files):
        print(BColors.cyan(f"--- Analyzing: {os.path.basename(sam_path)} ---"))
        with open(sam_path, 'r') as f_in:
            for line in f_in:
                if line.startswith('@'): continue
                
                try:
                    fields = line.strip().split('\t')
                    
                    if len(fields) < 10 or fields[9] == '*':
                        continue

                    rname = fields[2]
                    
                    if rname == '*':
                        continue

                    mutation_info = mutation_db.get(rname)
                    if not mutation_info:
                        continue
                    
                    if mutation_info['type'] not in allowed_variant_types:
                        continue

                    original_qname = fields[0]
                    suffixed_qname = f"{original_qname}_{file_idx + 1}"
                    
                    debug_row = None
                    if debug_mode:
                        debug_row = {
                            'query_id': suffixed_qname, 'ARO_matched': rname,
                            'status': 'N/A', 'reason': 'N/A', 'details': 'N/A', 'nucleotide_pid': 'N/A'
                        }

                    cigar = fields[5]
                    pos = int(fields[3])
                    ref_seq = db_sequences.get(rname)
                    
                    ref_span = get_reference_span_from_cigar(cigar)

                    if ref_span < min_aln_len:
                        if debug_mode:
                            debug_row['status'] = 'Rejected'
                            debug_row['reason'] = 'Alignment length below cutoff'
                            debug_row['details'] = f"Reference span {ref_span} < {min_aln_len}"
                            debug_writer.writerow(debug_row)
                        continue

                    if require_full_read_aln:
                        target_length = len(ref_seq) if ref_seq else 0
                        if not is_valid_alignment(cigar, pos, ref_span, target_length):
                            if debug_mode:
                                debug_row['status'] = 'Rejected'
                                debug_row['reason'] = 'Does not meet full read alignment criteria'
                                debug_row['details'] = f"CIGAR: {cigar}, POS: {pos}, REF_SPAN: {ref_span}, REF_LEN: {target_length}"
                                debug_writer.writerow(debug_row)
                            continue

                    optional_fields = fields[11:]
                    num_mismatches = get_nm_tag(optional_fields)
                    if num_mismatches is None: continue

                    nuc_pid = ((ref_span - num_mismatches) / ref_span) * 100.0 if ref_span > 0 else 0.0
                    if debug_mode:
                        debug_row['nucleotide_pid'] = f"{nuc_pid:.2f}"

                    query_seq = fields[9]

                    if ref_seq:
                        check_result = check_mutation_presence(cigar, pos, query_seq, ref_seq, mutation_info)
                        
                        if check_result['confirmed']:
                            confirmed_hits.append({
                                'query_id': suffixed_qname,
                                'ARO_matched': rname,
                                'mutation': mutation_info['full_mutation'],
                                'nucleotide_pid': f"{nuc_pid:.2f}",
                                'alignment_length': ref_span,
                                'position_on_ref': pos
                            })
                            alignment_str = format_fullseq_alignment_for_output(
                                cigar,
                                query_seq,
                                ref_seq,
                                mutation_info,
                                pos,
                                suffixed_qname,
                                rname
                            )
                            align_file.write(alignment_str)
                        
                        if debug_mode:
                            debug_row['status'] = 'Accepted' if check_result['confirmed'] else 'Rejected'
                            debug_row['reason'] = check_result['reason']
                            debug_row['details'] = check_result['details']
                            debug_writer.writerow(debug_row)

                except (IndexError, ValueError) as e:
                    print(BColors.yellow(f"Warning: Skipping malformed SAM line: {line.strip()} | Error: {e}"), file=sys.stderr)
                    continue
    
    if debug_file:
        debug_file.close()
    if align_file:
        align_file.close()

    print(BColors.green(f"\n--- Found a total of {len(confirmed_hits)} confirmed variant hits across all files. ---"))
    if not confirmed_hits:
        print(BColors.yellow("--- No variant hits were confirmed. Output file will contain only headers. ---"))

    try:
        with open(output_path, 'w', newline='') as f_out:
            header = ['query_id', 'ARO_matched', 'mutation', 'nucleotide_pid', 'alignment_length', 'position_on_ref']
            writer = csv.DictWriter(f_out, fieldnames=header, delimiter='\t')
            writer.writeheader()
            writer.writerows(confirmed_hits)
        print(BColors.green(f"--- Successfully wrote confirmed hits to: {output_path} ---"))
    except IOError as e:
        print(BColors.red(f"Error writing to output file {output_path}: {e}"), file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Analyzes SAM files to find and confirm specific resistance mutations."
    )
    parser.add_argument(
        "-i", "--input-files", 
        required=True, 
        help="Comma-delimited list of input SAM files."
    )
    parser.add_argument(
        "-d", "--db", 
        required=True, 
        help="Path to the reference database FASTA file, which contains mutation info in headers."
    )
    parser.add_argument(
        "--tmp-dir", 
        default=".", 
        help="Directory for output files. (Default: current directory)"
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix for the output TSV file (e.g., 'MyProject')."
    )
    parser.add_argument(
        "--min-aln-len",
        type=int,
        default=100,
        help="Minimum alignment length (in base pairs) on the reference for a hit to be reported. Default: 100"
    )
    parser.add_argument(
        "--require-full-read-aln",
        action="store_true",
        help="If specified, only count reads that align across their full length, allowing for clips only at the ends of a reference gene."
    )
    parser.add_argument(
        "--gene-types",
        default='V,R',
        help="Comma-delimited list of variant gene types to process (e.g., 'V,R,O'). Default: 'V,R'"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug mode, which creates a detailed TSV log file of all variant checks (`_variant_debug.tsv`)."
    )
    args = parser.parse_args()

    os.makedirs(args.tmp_dir, exist_ok=True)
    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    db_sequences = parse_fasta_db(args.db) 
    mutation_db = parse_mutation_info_from_db(db_sequences)

    allowed_variant_types = {t.strip().upper() for t in args.gene_types.split(',') if t.strip()}
    print(BColors.cyan(f"--- Processing variant types: {', '.join(sorted(list(allowed_variant_types)))} ---"))

    output_path = os.path.join(args.tmp_dir, f"{args.output_prefix}_variant_hits.tsv")
    debug_output_path = os.path.join(args.tmp_dir, f"{args.output_prefix}_variant_debug.tsv")
    alignment_output_path = os.path.join(args.tmp_dir, f"{args.output_prefix}_variant_alignments.txt")

    process_sam_files_for_variants(
        input_files_list, db_sequences, mutation_db, 
        output_path, args.debug, debug_output_path, alignment_output_path,
        args.min_aln_len, args.require_full_read_aln, allowed_variant_types
    )

    print(BColors.green("\n\n--- Variant Processing Complete ---"))

if __name__ == "__main__":
    main()