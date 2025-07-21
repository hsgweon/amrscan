#!/usr/bin/env python

import argparse
import os
import sys
import shutil
import re
import csv
from collections import defaultdict, Counter
# import io # Removed: Not needed without StreamRedirector

# --- BColors Class (consistent with other scripts) ---
class BColors:
    """A helper class to add color to terminal output."""
    HEADER = '\033[95m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

    # Detect if stdout is a tty, if not, disable colors
    enabled = sys.stdout.isatty()

    @classmethod
    def _colorize(cls, color_code, text):
        return f"{color_code}{text}{cls.ENDC}" if cls.enabled else text

    @classmethod
    def cyan(cls, text):
        return cls._colorize(cls.OKCYAN, text)

    @classmethod
    def green(cls, text):
        return cls._colorize(cls.OKGREEN, text)

    @classmethod
    def red(cls, text):
        # Errors may be printed to stderr, which might not be a TTY
        # Use stdout's tty status as the primary determinant
        return cls._colorize(cls.FAIL, text)

    @classmethod
    def yellow(cls, text):
        return cls._colorize(cls.WARNING, text)

# --- Global Constants (Python equivalent of C++ namespace) ---
class AMRScanFiles:
    COMMON_DB_PREFIX = "amrscan_DB_CARD"
    class Stage1:
        # LOG_SUFFIX = "" # Removed: No longer generating a log file
        README_SUFFIX = "README"
        FASTA_HOMOLOG = "homolog_model"
        FASTA_KNOCKOUT = "knockout_model"
        FASTA_VARIANT = "protein_variant_model"
        FASTA_OVEREXPRESSION = "protein_overexpression_model"
        FASTA_RRNA = "rRNA_gene_variant_model"
        FASTA_COMBINED = "all"
        METADATA = "metadata"
        METADATA_NON_VARIANT_IGNORED = "metadata_ignored_non_variant"
        METADATA_VARIANT_IGNORED = "metadata_ignored_variant"
    class Stage2:
        NR_METADATA = "NR_metadata"
        NR_FASTA = "NR_all"
        NR_FASTA_NUC = "NR_all_nuc"

# --- Data Structures (Python classes) ---
class Config:
    def __init__(self):
        self.path_snps_file = None
        self.path_protein_homolog_model = None
        self.path_protein_knockout_model = None
        self.path_protein_variant_model = None
        self.path_protein_overexpression_model = None
        self.path_rRNA_gene_variant_model = None
        self.path_aro_index = None
        self.output_dir_name = None
        self.card_version = None
        self.overwrite_outputs = False

class FastaEntry:
    def __init__(self, header="", sequence=""):
        self.header = header
        self.sequence = sequence

class ProteinFastaHeaderInfo:
    def __init__(self):
        self.aro_number = "N/A"
        self.name = "N/A"
        self.card_short_name = "N/A"
        self.parsed_successfully = False

class SNPSFileLine:
    def __init__(self):
        self.accession_aro = ""
        self.name = ""
        self.model_type = ""
        self.parameter_type = ""
        self.mutations_str = ""
        self.card_short_name = ""
        self.source = ""
        self.citation = ""
        self.original_line = ""
        self.line_number = 0
        self.is_protein_pathway_candidate = False
        self.is_rrna_pathway_candidate = False

class ParsedProteinMutation:
    def __init__(self):
        self.original_aa = ''
        self.new_aa = ''
        self.position = 0
        self.is_frameshift = False
        self.is_nonsense = False
        self.original_token = ''

class ParsedNucleotideMutation:
    def __init__(self):
        self.original_base = ''
        self.new_base = ''
        self.position = 0
        self.original_token = ''

class MetadataEntry:
    def __init__(self, line=None, is_header=False):
        self.sequence_id = ""
        self.aro_number = ""
        self.card_short_name = ""
        self.model_type = ""
        self.parameter_type = ""
        self.mutation_string_in_log = ""
        self.nucleotide_mutation_position_s = ""
        self.seq_nuc_length = 0
        self.original_line = line if line is not None else ""

        if line and not is_header:
            cols = split_string(line, '\t')
            if len(cols) >= 8:
                self.sequence_id = cols[0]
                self.aro_number = cols[1]
                self.card_short_name = cols[2]
                self.model_type = cols[3]
                self.parameter_type = cols[4]
                self.mutation_string_in_log = cols[5]
                self.nucleotide_mutation_position_s = cols[6]
                try: self.seq_nuc_length = int(cols[7])
                except (ValueError, IndexError): self.seq_nuc_length = 0
            else:
                self.sequence_id = "ERROR_PARSING"

class AroIndexEntry:
    def __init__(self, line=None, line_num=0):
        self.aro_accession = "N/A"
        self.cvterm_id = ""
        self.model_sequence_id = ""
        self.model_id = ""
        self.model_name = ""
        self.aro_name = ""
        self.protein_accession = ""
        self.dna_accession = ""
        self.amr_gene_family = ""
        self.drug_class = ""
        self.resistance_mechanism = ""
        self.card_short_name = ""

        if line:
            cols = split_string(line, '\t')
            if len(cols) >= 12:
                colon_pos = cols[0].rfind(':')
                self.aro_accession = cols[0][colon_pos + 1:] if colon_pos != -1 and colon_pos + 1 < len(cols[0]) else "N/A"
                self.cvterm_id = cols[1]
                self.model_sequence_id = cols[2]
                self.model_id = cols[3]
                self.model_name = cols[4]
                self.aro_name = cols[5]
                self.protein_accession = cols[6]
                self.dna_accession = cols[7]
                self.amr_gene_family = cols[8]
                self.drug_class = cols[9]
                self.resistance_mechanism = cols[10]
                self.card_short_name = cols[11]

# --- StreamRedirector Class (Removed) ---

# --- Helper Functions ---
def trim_string(s):
    return s.strip()

def split_string(s, delimiter):
    return [trim_string(token) for token in s.split(delimiter)]

def to_uppercase_str(s):
    return s.upper()

def file_exists(path):
    return os.path.exists(path) and os.path.isfile(path)

def format_output_filename(suffix_part, version, ext):
    # LOG_SUFFIX is no longer relevant here as log file is removed
    if not suffix_part: # This condition was for LOG_SUFFIX, now it means no suffix
        return f"{AMRScanFiles.COMMON_DB_PREFIX}_v{version}{ext}"
    else:
        return f"{AMRScanFiles.COMMON_DB_PREFIX}_v{version}_{suffix_part}{ext}"

def sanitize_for_id(s):
    if not s or s == "NA" or s == "Wildtype": # Treat "NA" and "Wildtype" as empty for ID generation
        return ""
    s = re.sub(r'[ ()\[\]/|*,.;:]', '_', s)
    s = re.sub(r'__+', '_', s) # Replace multiple underscores with single
    s = s.strip('_') # Remove leading/trailing underscores
    return s

def get_model_code(model_type):
    if "homolog" in model_type: return 'H'
    if "knockout" in model_type: return 'K'
    if "protein variant" in model_type: return 'V'
    if "overexpression" in model_type: return 'O'
    if "rRNA" in model_type: return 'R'
    return 'U'

def generate_unique_id(aro, name, model_type, mutation, id_tracker):
    mut_sanitized = sanitize_for_id(mutation) # sanitize_for_id now handles "NA" and "Wildtype" to empty string
    
    # Conditionally add mutation part only if it's not empty after sanitization
    mutation_part = f"__{mut_sanitized}" if mut_sanitized else ""

    base_id = f"ARO_{sanitize_for_id(aro)}__{sanitize_for_id(name)}__{get_model_code(model_type)}{mutation_part}"

    count = id_tracker[base_id]
    id_tracker[base_id] += 1
    return f"{base_id}_{count}" if count > 0 else base_id

# --- Codon Translation (consistent with other scripts) ---
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Returns a representative codon for each amino acid. Not exhaustive.
REVERSE_CODON_TABLE = {
    'A': "GCT", 'C': "TGT", 'D': "GAT", 'E': "GAA", 'F': "TTT",
    'G': "GGT", 'H': "CAT", 'I': "ATT", 'K': "AAA", 'L': "TTA",
    'M': "ATG", 'N': "AAT", 'P': "CCT", 'Q': "CAA", 'R': "CGT",
    'S': "TCT", 'T': "ACT", 'V': "GTT", 'W': "TGG", 'Y': "TAT",
    '*': "TAA"
}

def translate_codon(codon):
    """Translates a single codon. Handles gaps and ambiguous bases."""
    if len(codon) != 3 or '-' in codon or 'N' in codon.upper():
        return 'X'
    return CODON_TABLE.get(codon.upper(), 'X')

def translate_nucleotide_sequence(nuc_seq):
    """
    Translates a full nucleotide sequence to protein.
    Handles 'U' to 'T' conversion, skips non-ACGT characters, and stops at stop codons.
    Returns (protein_sequence, error_message).
    """
    prot_seq = []
    err_msgs = []
    valid_nuc_seq = ""

    for char in to_uppercase_str(nuc_seq):
        if char in "ACGT":
            valid_nuc_seq += char
        elif char == 'U':
            valid_nuc_seq += 'T'
        elif not char.isspace():
            err_msgs.append(f"Invalid char '{char}'.")

    if len(valid_nuc_seq) < 3:
        err_msgs.append("Sequence too short for translation.")
        return "", "; ".join(err_msgs)

    for i in range(0, len(valid_nuc_seq) - 2, 3):
        codon = valid_nuc_seq[i:i+3]
        aa = translate_codon(codon)
        prot_seq.append(aa)
        if aa == '*':
            break # Stop translation at the first stop codon

    return "".join(prot_seq), "; ".join(err_msgs)

def join_set_to_string(s_set, delimiter):
    if not s_set:
        return "NA" # Keep "NA" for metadata fields if no value
    return delimiter.join(sorted(list(s_set)))

def join_vector_to_delimited_string(s_vec, delimiter):
    if not s_vec:
        return "NA" # Keep "NA" for metadata fields if no value
    return delimiter.join(s_vec)

def read_config(config_path):
    config = Config()
    config_map = {}
    try:
        with open(config_path, 'r') as f:
            for line in f:
                line = trim_string(line)
                if not line or line.startswith('#'):
                    continue
                if ':' in line:
                    key, value = line.split(':', 1)
                    config_map[trim_string(key)] = trim_string(value)
    except IOError as e:
        raise RuntimeError(f"Error: Could not open configuration file: {config_path} - {e}")

    def get_config_value(key):
        if key not in config_map or not config_map[key]:
            raise RuntimeError(f"Mandatory config key missing or empty: {key}")
        return config_map[key]

    config.path_snps_file = get_config_value("PATH_snps_file")
    config.path_protein_homolog_model = get_config_value("PATH_protein_homolog_model")
    config.path_protein_knockout_model = get_config_value("PATH_protein_knockout_model")
    config.path_protein_variant_model = get_config_value("PATH_protein_variant_model")
    config.path_protein_overexpression_model = get_config_value("PATH_protein_overexpression_model")
    config.path_rRNA_gene_variant_model = get_config_value("PATH_rRNA_gene_variant_model")
    config.path_aro_index = get_config_value("PATH_aro_index")
    config.card_version = get_config_value("CARD_VERSION")

    return config

def display_help():
    print(BColors.cyan("AMRScan_DB_CARD_prep - Processes CARD data into a non-redundant, annotated database.\n"))
    print(BColors.cyan("Usage: AMRScan_DB_CARD_prep -c <config> -d <out_dir> [--overwrite]\n"))
    print(BColors.cyan("Options:"))
    print(BColors.cyan("  -c, --config <path>      Path to the configuration file (mandatory)."))
    print(BColors.cyan("  -d, --output-dir <name>  Name of the output directory (mandatory)."))
    print(BColors.cyan("  --overwrite              Delete and recreate output directory if it exists."))
    print(BColors.cyan("  -h, --help               Display this help message."))

def setup_output_directory(dir_path, overwrite):
    if os.path.exists(dir_path):
        if not os.path.isdir(dir_path):
            raise RuntimeError(f"Error: Output path exists but is not a directory: {dir_path}")
        if overwrite:
            print(BColors.yellow("Info: Removing existing output directory."))
            shutil.rmtree(dir_path)
        else:
            raise RuntimeError(f"Error: Output directory '{dir_path}' exists. Use --overwrite to replace it.")
    os.makedirs(dir_path)
    print(BColors.green(f"Info: Created output directory: {dir_path}"))
    return True

def read_next_fasta_entry(fasta_file_handle):
    entry = FastaEntry()
    line = fasta_file_handle.readline()
    while line:
        line = trim_string(line)
        if line and line.startswith('>'):
            entry.header = line
            break
        line = fasta_file_handle.readline()

    if not entry.header:
        return None

    while True:
        pos = fasta_file_handle.tell() # Remember current position
        line = fasta_file_handle.readline()
        if not line or line.startswith('>'):
            fasta_file_handle.seek(pos) # Go back to start of next entry or EOF
            break
        entry.sequence += trim_string(line)

    entry.sequence = to_uppercase_str(entry.sequence)
    return entry

def parse_protein_fasta_header(header_line):
    info = ProteinFastaHeaderInfo()
    content = header_line[1:] # Remove '>'
    
    match_aro = re.search(r'\|ARO:(\d+)', content)
    if match_aro:
        info.aro_number = match_aro.group(1)

    last_pipe_pos = content.rfind('|')
    name_part = trim_string(content[last_pipe_pos + 1:]) if last_pipe_pos != -1 else content
    
    if name_part:
        info.name = name_part
        # Find end of short name (space, [, or end of string)
        end_pos = len(name_part)
        space_pos = name_part.find(' ')
        bracket_pos = name_part.find('[')
        
        if space_pos != -1:
            end_pos = min(end_pos, space_pos)
        if bracket_pos != -1:
            end_pos = min(end_pos, bracket_pos)
            
        info.card_short_name = trim_string(name_part[:end_pos])

    info.parsed_successfully = (info.aro_number != "N/A" or info.card_short_name != "N/A")
    return info

def process_protein_model_fasta(in_path, suffix, model_type, metadata_file, ignored_file, db_dir, version, id_tracker, sequence_to_nucleotide_map):
    try:
        in_fasta = open(in_path, 'r')
    except IOError:
        ignored_file.write(f"SYSTEM_ERROR\tCould not open: {in_path}\n")
        return

    out_prot_path = os.path.join(db_dir, format_output_filename(suffix, version, ".fasta"))
    out_nuc_path = os.path.join(db_dir, format_output_filename(suffix + "_nuc", version, ".fasta"))

    try:
        out_prot_fasta = open(out_prot_path, 'w')
        out_nuc_fasta = open(out_nuc_path, 'w')
    except IOError as e:
        ignored_file.write(f"SYSTEM_ERROR\tCould not create output files for {suffix}: {e}\n")
        in_fasta.close()
        return

    while True:
        entry = read_next_fasta_entry(in_fasta)
        if not entry:
            break

        prot_seq, err = translate_nucleotide_sequence(entry.sequence)
        h_info = parse_protein_fasta_header(entry.header)

        if h_info.parsed_successfully and prot_seq and not err:
            # Pass empty string for mutation, so generate_unique_id omits the mutation part
            seq_id = generate_unique_id(h_info.aro_number, h_info.card_short_name, model_type, "", id_tracker)
            out_prot_fasta.write(f">{seq_id}\n{prot_seq}\n")
            out_nuc_fasta.write(f">{seq_id}\n{entry.sequence}\n")
            metadata_file.write(f"{seq_id}\t{h_info.aro_number}\t{h_info.card_short_name}\t"
                                f"{model_type}\tNA\tNA\tNA\t{len(entry.sequence)}\n")
            sequence_to_nucleotide_map[prot_seq] = entry.sequence
        else:
            ignored_file.write(f"{entry.header}\tParse/translation failed: {err}\n")
    
    in_fasta.close()
    out_prot_fasta.close()
    out_nuc_fasta.close()

def load_fasta_to_map(path):
    fasta_map = {}
    try:
        with open(path, 'r') as f:
            while True:
                entry = read_next_fasta_entry(f)
                if not entry:
                    break
                match_aro = re.search(r'\|ARO:(\d+)', entry.header)
                if match_aro:
                    fasta_map[match_aro.group(1)] = entry
    except IOError:
        return False, fasta_map
    return True, fasta_map

def parse_snps_file_line(line, line_num):
    entry = SNPSFileLine()
    entry.original_line = line
    entry.line_number = line_num
    cols = split_string(line, '\t')
    if len(cols) >= 6:
        entry.accession_aro = cols[0]
        entry.name = cols[1]
        entry.model_type = cols[2]
        entry.parameter_type = cols[3]
        entry.mutations_str = cols[4]
        entry.card_short_name = cols[5]
        
        is_prot = "protein" in entry.model_type
        is_rrna = "rRNA" in entry.model_type
        is_var = "variant" in entry.parameter_type or "mutation" in entry.parameter_type
        
        entry.is_protein_pathway_candidate = is_prot and is_var
        entry.is_rrna_pathway_candidate = is_rrna and is_var
    return entry

def parse_protein_mutations_string(mut_str):
    muts = []
    if not mut_str or trim_string(mut_str) == "NA":
        return True, muts # Empty or NA string is valid (wildtype)

    for token in split_string(mut_str, ','):
        if not token: continue
        p_mut = ParsedProteinMutation()
    #     p_mut.original_token = token # This line was causing an error in the C++ version, but not in Python. Keeping it for consistency.
        p_mut.original_token = token # Re-added for consistency with C++ struct

        # Frameshift: (AA)(Pos)fs
        match_fs = re.match(r'([A-Z])(\d+)fs', token)
        if match_fs:
            p_mut.original_aa = match_fs.group(1)[0]
            p_mut.position = int(match_fs.group(2))
            p_mut.is_frameshift = True
            muts.append(p_mut)
            continue

        # Nonsense: (AA)(Pos)Ter
        match_ter = re.match(r'([A-Z])(\d+)Ter', token)
        if match_ter:
            p_mut.original_aa = match_ter.group(1)[0]
            p_mut.position = int(match_ter.group(2))
            p_mut.new_aa = '*'
            p_mut.is_nonsense = True
            muts.append(p_mut)
            continue

        # Substitution: (OriginalAA)(Pos)(NewAA)
        match_sub = re.match(r'([A-Z])(\d+)([A-Z*])', token)
        if match_sub:
            p_mut.original_aa = match_sub.group(1)[0]
            p_mut.position = int(match_sub.group(2))
            p_mut.new_aa = match_sub.group(3)[0]
            muts.append(p_mut)
            continue
        
        return False, [] # Failed to parse token

    return True, muts

def parse_nucleotide_mutations_string(mut_str):
    muts = []
    if not mut_str or trim_string(mut_str) == "NA":
        return True, muts # Empty or NA string is valid (wildtype)

    for token in split_string(mut_str, ','):
        if not token: continue
        n_mut = ParsedNucleotideMutation()
        n_mut.original_token = token
        
        # Nucleotide substitution: (OriginalBase)(Pos)(NewBase)
        match_nuc = re.match(r'([acgtuACGTU])(\d+)([acgtuACGTU])', token)
        if match_nuc:
            n_mut.original_base = to_uppercase_str(match_nuc.group(1))[0]
            n_mut.position = int(match_nuc.group(2))
            n_mut.new_base = to_uppercase_str(match_nuc.group(3))[0]
            if n_mut.original_base == 'U': n_mut.original_base = 'T'
            if n_mut.new_base == 'U': n_mut.new_base = 'T'
            muts.append(n_mut)
            continue
        
        return False, [] # Failed to parse token

    return True, muts

def process_protein_variant_entry(snps, source_map, out_prot_fasta, out_nuc_fasta, metadata_file, ignored_file, id_tracker, sequence_to_nucleotide_map):
    source_entry = source_map.get(snps.accession_aro)
    if not source_entry or not source_entry.sequence:
        ignored_file.write(f"{snps.original_line}\tError: ARO not found or empty sequence in source FASTA\n")
        return

    def write_entry(prot_seq, nuc_seq, mut_tag, nuc_pos):
        if not prot_seq: return
        seq_id = generate_unique_id(snps.accession_aro, snps.card_short_name, snps.model_type, mut_tag, id_tracker)
        out_prot_fasta.write(f">{seq_id}\n{prot_seq}\n")
        out_nuc_fasta.write(f">{seq_id}\n{nuc_seq}\n")
        metadata_file.write(f"{seq_id}\t{snps.accession_aro}\t{snps.card_short_name}\t"
                            f"{snps.model_type}\t{snps.parameter_type}\t{mut_tag}\t"
                            f"{nuc_pos}\t{len(nuc_seq)}\n")
        sequence_to_nucleotide_map[prot_seq] = nuc_seq

    parse_success, mutations = parse_protein_mutations_string(snps.mutations_str)
    if not parse_success:
        ignored_file.write(f"{snps.original_line}\tError: Could not parse protein mutations string\n")
        return

    if not mutations: # Wildtype
        prot_seq, err = translate_nucleotide_sequence(source_entry.sequence)
        if not prot_seq or err:
            ignored_file.write(f"{snps.original_line}\tError: Wildtype translation failed: {err}\n")
            return
        # Pass empty string for mutation, so generate_unique_id omits the mutation part
        write_entry(prot_seq, source_entry.sequence, "", "NA")
        return

    # Handle frameshift mutations
    if mutations[0].is_frameshift:
        fs = mutations[0]
        nuc_pos_0 = (fs.position - 1) * 3
        for del_len in [1, 2]: # Simulate 1bp and 2bp deletions for frameshift
            if nuc_pos_0 >= 0 and nuc_pos_0 + del_len <= len(source_entry.sequence):
                nt_del = list(source_entry.sequence) # Convert to list for mutable string
                del nt_del[nuc_pos_0 : nuc_pos_0 + del_len]
                nt_del_str = "".join(nt_del)
                
                prot_del, fs_err = translate_nucleotide_sequence(nt_del_str)
                if prot_del and not fs_err:
                    write_entry(prot_del, nt_del_str, f"{fs.original_token}_{del_len}del", str(nuc_pos_0 + 1))
                else:
                    ignored_file.write(f"{snps.original_line}\tError: Frameshift translation failed for {fs.original_token}_{del_len}del: {fs_err}\n")
            else:
                ignored_file.write(f"{snps.original_line}\tError: Frameshift position out of bounds for {fs.original_token}\n")
        return

    # Handle substitution/nonsense mutations
    current_nuc_list = list(source_entry.sequence)
    nuc_pos_vec = []
    
    orig_prot, err_trans = translate_nucleotide_sequence("".join(current_nuc_list))
    if not orig_prot:
        ignored_file.write(f"{snps.original_line}\tError: Initial translation failed for substitution: {err_trans}\n")
        return

    for mut in mutations:
        codon_start_pos = (mut.position - 1) * 3
        if mut.position <= 0 or codon_start_pos + 3 > len(current_nuc_list):
            ignored_file.write(f"{snps.original_line}\tError: Mutation position out of bounds for {mut.original_token}\n")
            return

        # Verify original AA at position
        current_prot_check, _ = translate_nucleotide_sequence("".join(current_nuc_list))
        if mut.position > len(current_prot_check) or current_prot_check[mut.position-1] != mut.original_aa:
            # Allow if the current AA is already the new AA (e.g., multiple mutations on same codon)
            if mut.position > len(current_prot_check) or current_prot_check[mut.position-1] != mut.new_aa:
                ignored_file.write(f"{snps.original_line}\tError: AA mismatch for {mut.original_token}. Expected '{mut.original_aa}', found '{current_prot_check[mut.position-1] if mut.position <= len(current_prot_check) else 'N/A'}'\n")
                return

        if mut.new_aa in REVERSE_CODON_TABLE:
            new_codon = REVERSE_CODON_TABLE[mut.new_aa]
            for i in range(3):
                current_nuc_list[codon_start_pos + i] = new_codon[i]
            nuc_pos_vec.append(str(codon_start_pos + 1))
        else:
            ignored_file.write(f"{snps.original_line}\tError: No reverse codon found for new AA '{mut.new_aa}' in {mut.original_token}\n")
            return

    final_nuc_seq = "".join(current_nuc_list)
    final_prot, final_err = translate_nucleotide_sequence(final_nuc_seq)
    if not final_prot or final_err:
        ignored_file.write(f"{snps.original_line}\tError: Final translation failed for substitution: {final_err}\n")
        return
    
    write_entry(final_prot, final_nuc_seq, snps.mutations_str, join_vector_to_delimited_string(nuc_pos_vec, ","))

def process_rrna_variant_entry(snps, rrna_map, out_fasta, metadata_file, ignored_file, id_tracker, sequence_to_nucleotide_map):
    source_entry = rrna_map.get(snps.accession_aro)
    if not source_entry or not source_entry.sequence:
        ignored_file.write(f"{snps.original_line}\tError: ARO not found or empty sequence in source FASTA\n")
        return

    def write_entry(nuc_seq, mut_tag, nuc_pos):
        seq_id = generate_unique_id(snps.accession_aro, snps.card_short_name, snps.model_type, mut_tag, id_tracker)
        out_fasta.write(f">{seq_id}\n{nuc_seq}\n")
        metadata_file.write(f"{seq_id}\t{snps.accession_aro}\t{snps.card_short_name}\t"
                            f"{snps.model_type}\t{snps.parameter_type}\t{mut_tag}\t"
                            f"{nuc_pos}\t{len(nuc_seq)}\n")
        sequence_to_nucleotide_map[nuc_seq] = nuc_seq # For rRNA, nuc_seq is the "protein" equivalent

    parse_success, mutations = parse_nucleotide_mutations_string(snps.mutations_str)
    if not parse_success:
        ignored_file.write(f"{snps.original_line}\tError: Could not parse rRNA mutations string\n")
        return

    if not mutations: # Wildtype
        # Pass empty string for mutation, so generate_unique_id omits the mutation part
        write_entry(source_entry.sequence, "", "NA")
        return

    current_nuc_list = list(source_entry.sequence)
    nuc_pos_vec = []
    for mut in mutations:
        if mut.position <= 0 or mut.position > len(current_nuc_list):
            ignored_file.write(f"{snps.original_line}\tError: rRNA Mutation position out of bounds for {mut.original_token}\n")
            return
        
        # Check if original base matches, or if it's already the new base (multiple mutations on same pos)
        if current_nuc_list[mut.position - 1] != mut.original_base and current_nuc_list[mut.position - 1] != mut.new_base:
            ignored_file.write(f"{snps.original_line}\tError: rRNA Mutation mismatch for {mut.original_token}. Expected '{mut.original_base}', found '{current_nuc_list[mut.position-1]}'\n")
            return
        
        current_nuc_list[mut.position - 1] = mut.new_base
        nuc_pos_vec.append(str(mut.position))
    
    final_nuc_seq = "".join(current_nuc_list)
    write_entry(final_nuc_seq, snps.mutations_str, join_vector_to_delimited_string(nuc_pos_vec, ","))

def load_all_metadata(path):
    entries = []
    try:
        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader) # Skip header
            for row in reader:
                entries.append(MetadataEntry('\t'.join(row)))
    except IOError:
        raise RuntimeError(f"Failed to load Stage 1 metadata from {path}")
    return entries

def load_all_fasta_sequences(path):
    fasta_map = {}
    try:
        with open(path, 'r') as f:
            while True:
                entry = read_next_fasta_entry(f)
                if not entry:
                    break
                if entry.header and len(entry.header) > 1:
                    fasta_map[trim_string(entry.header[1:])] = entry
    except IOError:
        raise RuntimeError(f"Failed to load combined Stage 1 FASTA from {path}")
    return fasta_map

def load_aro_index_file(path):
    aro_index_map = {}
    try:
        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader) # Skip header
            for i, row in enumerate(reader):
                entry = AroIndexEntry('\t'.join(row), i + 1)
                if entry.aro_accession != "N/A":
                    aro_index_map[entry.aro_accession] = entry
    except IOError:
        raise RuntimeError(f"Failed to load ARO index file from {path}")
    return aro_index_map

def run_first_stage(config, db_dir, sequence_to_nucleotide_map):
    print(BColors.cyan("\n--- Stage 1: Processing CARD Models and SNPs ---"))
    
    stage1_id_tracker = defaultdict(int)

    def open_output_stream(suffix, ext, header):
        path = os.path.join(db_dir, format_output_filename(suffix, config.card_version, ext))
        stream = open(path, 'w')
        if header:
            stream.write(header + '\n')
        return stream

    metadata_file = open_output_stream(AMRScanFiles.Stage1.METADATA, ".txt", "Sequence_ID\tARO_Number\tCard_Short_Name\tModel_Type\tParameter_Type\tMutation_String\tNucleotide_Mutation_Position\tSeqNucLength")
    non_variant_ignored = open_output_stream(AMRScanFiles.Stage1.METADATA_NON_VARIANT_IGNORED, ".txt", "Original_Header\tIssue_Description")
    variant_ignored = open_output_stream(AMRScanFiles.Stage1.METADATA_VARIANT_IGNORED, ".txt", "Original_snps.txt_Line\tIssue_Description")

    print(BColors.cyan("Processing homolog models..."))
    process_protein_model_fasta(config.path_protein_homolog_model, AMRScanFiles.Stage1.FASTA_HOMOLOG, "homolog", metadata_file, non_variant_ignored, db_dir, config.card_version, stage1_id_tracker, sequence_to_nucleotide_map)
    print(BColors.green("Homolog models processed."))

    print(BColors.cyan("Processing knockout models..."))
    process_protein_model_fasta(config.path_protein_knockout_model, AMRScanFiles.Stage1.FASTA_KNOCKOUT, "knockout", metadata_file, non_variant_ignored, db_dir, config.card_version, stage1_id_tracker, sequence_to_nucleotide_map)
    print(BColors.green("Knockout models processed."))

    print(BColors.cyan("Loading variant FASTA files..."))
    success, protein_variant_map = load_fasta_to_map(config.path_protein_variant_model)
    if not success: raise RuntimeError("Failed to load Protein Variant FASTA")
    success, protein_overexpression_map = load_fasta_to_map(config.path_protein_overexpression_model)
    if not success: raise RuntimeError("Failed to load Protein Overexpression FASTA")
    success, rrna_variant_map = load_fasta_to_map(config.path_rRNA_gene_variant_model)
    if not success: raise RuntimeError("Failed to load rRNA Variant FASTA")
    print(BColors.green("Variant FASTA files loaded."))

    prot_var_fasta = open_output_stream(AMRScanFiles.Stage1.FASTA_VARIANT, ".fasta", "")
    prot_var_nuc_fasta = open_output_stream(AMRScanFiles.Stage1.FASTA_VARIANT + "_nuc", ".fasta", "")
    prot_over_fasta = open_output_stream(AMRScanFiles.Stage1.FASTA_OVEREXPRESSION, ".fasta", "")
    prot_over_nuc_fasta = open_output_stream(AMRScanFiles.Stage1.FASTA_OVEREXPRESSION + "_nuc", ".fasta", "")
    rrna_fasta = open_output_stream(AMRScanFiles.Stage1.FASTA_RRNA, ".fasta", "")

    print(BColors.cyan("Processing snps.txt for variants..."))
    try:
        with open(config.path_snps_file, 'r') as snps_file:
            csv_reader = csv.reader(snps_file, delimiter='\t')
            next(csv_reader) # Skip header
            for i, row in enumerate(csv_reader):
                line = '\t'.join(row)
                entry = parse_snps_file_line(line, i + 1)
                if entry.is_protein_pathway_candidate:
                    if "variant" in entry.model_type:
                        process_protein_variant_entry(entry, protein_variant_map, prot_var_fasta, prot_var_nuc_fasta, metadata_file, variant_ignored, stage1_id_tracker, sequence_to_nucleotide_map)
                    elif "overexpression" in entry.model_type:
                        process_protein_variant_entry(entry, protein_overexpression_map, prot_over_fasta, prot_over_nuc_fasta, metadata_file, variant_ignored, stage1_id_tracker, sequence_to_nucleotide_map)
                elif entry.is_rrna_pathway_candidate:
                    process_rrna_variant_entry(entry, rrna_variant_map, rrna_fasta, metadata_file, variant_ignored, stage1_id_tracker, sequence_to_nucleotide_map)
                elif entry.accession_aro != "": # If it's a valid line but not a candidate
                    variant_ignored.write(f"{entry.original_line}\tSkipped: Model/Parameter Type not applicable.\n")
    except IOError:
        raise RuntimeError(f"Failed to open snps file: {config.path_snps_file}")
    print(BColors.green("snps.txt processed."))

    metadata_file.close()
    non_variant_ignored.close()
    variant_ignored.close()
    prot_var_fasta.close()
    prot_var_nuc_fasta.close()
    prot_over_fasta.close()
    prot_over_nuc_fasta.close()
    rrna_fasta.close()
    print(BColors.green("--- Stage 1 Complete ---"))

def run_second_stage(config, db_dir, sequence_to_nucleotide_map):
    print(BColors.cyan("\n--- Stage 2: Removing Redundancy and Finalizing Outputs ---"))

    combined_fasta_path = os.path.join(db_dir, format_output_filename(AMRScanFiles.Stage1.FASTA_COMBINED, config.card_version, ".fasta"))
    
    # Combine all individual FASTA files
    with open(combined_fasta_path, 'w') as combined_stream:
        for suffix in [AMRScanFiles.Stage1.FASTA_HOMOLOG, AMRScanFiles.Stage1.FASTA_KNOCKOUT, AMRScanFiles.Stage1.FASTA_VARIANT, AMRScanFiles.Stage1.FASTA_OVEREXPRESSION, AMRScanFiles.Stage1.FASTA_RRNA]:
            individual_fasta_path = os.path.join(db_dir, format_output_filename(suffix, config.card_version, ".fasta"))
            if os.path.exists(individual_fasta_path):
                with open(individual_fasta_path, 'r') as individual_stream:
                    shutil.copyfileobj(individual_stream, combined_stream)
    print(BColors.green("Combined all Stage 1 FASTA files."))

    all_metadata = load_all_metadata(os.path.join(db_dir, format_output_filename(AMRScanFiles.Stage1.METADATA, config.card_version, ".txt")))
    id_to_fasta = load_all_fasta_sequences(combined_fasta_path)
    aro_index_map = load_aro_index_file(config.path_aro_index)

    id_to_metadata = {}
    sequence_to_ids = defaultdict(list)
    for meta in all_metadata:
        if meta.sequence_id == "ERROR_PARSING":
            continue
        id_to_metadata[meta.sequence_id] = meta
        if meta.sequence_id in id_to_fasta:
            sequence_to_ids[id_to_fasta[meta.sequence_id].sequence].append(meta.sequence_id)
    print(BColors.green(f"Loaded {len(id_to_metadata)} metadata entries and {len(id_to_fasta)} FASTA entries."))

    nr_meta_path = os.path.join(db_dir, format_output_filename(AMRScanFiles.Stage2.NR_METADATA, config.card_version, ".txt"))
    nr_fasta_path = os.path.join(db_dir, format_output_filename(AMRScanFiles.Stage2.NR_FASTA, config.card_version, ".fasta"))
    nr_nuc_fasta_path = os.path.join(db_dir, format_output_filename(AMRScanFiles.Stage2.NR_FASTA_NUC, config.card_version, ".fasta"))

    with open(nr_meta_path, 'w', newline='') as nr_meta_stream, \
         open(nr_fasta_path, 'w') as nr_fasta_stream, \
         open(nr_nuc_fasta_path, 'w') as nr_nuc_fasta_stream:
        
        nr_meta_writer = csv.writer(nr_meta_stream, delimiter='\t')
        nr_meta_writer.writerow(["Sequence_ID", "ARO_Number", "Card_Short_Name", "Model_Type", "Parameter_Type",
                                 "Mutation_String", "Nucleotide_Mutation_Position", "SeqNucLength",
                                 "Consolidated_Stage1_IDs", "Consolidated_ARO_Numbers", "Consolidated_Mutation_Strings",
                                 "Consolidated_Parameter_Types", "Redundancy_Count", "AMR_Gene_Family", "ARO_Name",
                                 "Drug_Class", "Resistance_Mechanism"])
        
        final_id_tracker = defaultdict(int)

        for unique_sequence, ids in sequence_to_ids.items():
            if not ids: continue

            rep_meta = id_to_metadata[ids[0]] # Representative metadata from the first ID

            consolidated_aros = set()
            consolidated_mutations = set()
            consolidated_params = set()
            consolidated_families = set()
            consolidated_aro_names = set()
            consolidated_drugs = set()
            consolidated_mechs = set()

            for id_val in ids:
                meta = id_to_metadata[id_val]
                if meta.aro_number != "NA":
                    consolidated_aros.add(meta.aro_number)
                    if meta.aro_number in aro_index_map:
                        aro_data = aro_index_map[meta.aro_number]
                        if aro_data.amr_gene_family and aro_data.amr_gene_family != "NA": consolidated_families.add(aro_data.amr_gene_family)
                        if aro_data.aro_name and aro_data.aro_name != "NA": consolidated_aro_names.add(aro_data.aro_name)
                        if aro_data.drug_class and aro_data.drug_class != "NA": consolidated_drugs.add(aro_data.drug_class)
                        if aro_data.resistance_mechanism and aro_data.resistance_mechanism != "NA": consolidated_mechs.add(aro_data.resistance_mechanism)
                if meta.mutation_string_in_log not in ["NA", "Wildtype"]:
                    consolidated_mutations.add(meta.mutation_string_in_log)
                if meta.parameter_type != "NA":
                    consolidated_params.add(meta.parameter_type)

            consolidated_aros_str = join_set_to_string(consolidated_aros, ";")
            consolidated_muts_str = join_set_to_string(consolidated_mutations, ";")

            # Sanitize consolidated mutations for header, will be empty if "NA" or "Wildtype"
            consolidated_muts_sanitized_for_header = sanitize_for_id(consolidated_muts_str)
            mutation_part_for_header = f"__{consolidated_muts_sanitized_for_header}" if consolidated_muts_sanitized_for_header else ""

            base_id = (f"ARO_{sanitize_for_id(consolidated_aros_str)}__{sanitize_for_id(rep_meta.card_short_name)}__"
                       f"{get_model_code(rep_meta.model_type)}{mutation_part_for_header}")
            
            count = final_id_tracker[base_id]
            final_id_tracker[base_id] += 1
            final_id = f"{base_id}_{count}" if count > 0 else base_id

            if count > 0:
                sys.stderr.write(BColors.yellow(f"Warning: Generated identical base ID '{base_id}' for a new unique sequence.\n"
                                                f"         This usually means different mutations resulted in the same final sequence.\n"
                                                f"         Appending suffix to ensure unique ID: {final_id}\n"))

            nr_fasta_stream.write(f">{final_id}\n{unique_sequence}\n")
            
            nuc_seq_for_protein = sequence_to_nucleotide_map.get(unique_sequence)
            if nuc_seq_for_protein:
                nr_nuc_fasta_stream.write(f">{final_id}\n{nuc_seq_for_protein}\n")
            else:
                sys.stderr.write(BColors.yellow(f"Warning: Could not find corresponding nucleotide sequence for final ID: {final_id}\n"))
            
            nr_meta_writer.writerow([
                final_id, rep_meta.aro_number, rep_meta.card_short_name, rep_meta.model_type, rep_meta.parameter_type,
                rep_meta.mutation_string_in_log, rep_meta.nucleotide_mutation_position_s, rep_meta.seq_nuc_length,
                join_vector_to_delimited_string(ids, ";"),
                consolidated_aros_str,
                (consolidated_muts_str if consolidated_muts_str != "NA" else ""), # Empty string if no mutations
                join_set_to_string(consolidated_params, ";"),
                len(ids),
                join_set_to_string(consolidated_families, ";"),
                join_set_to_string(consolidated_aro_names, ";"),
                join_set_to_string(consolidated_drugs, ";"),
                join_set_to_string(consolidated_mechs, ";")
            ])
    print(BColors.green("--- Stage 2 Complete ---"))

def write_readme_file(config, db_dir):
    readme_path = os.path.join(db_dir, format_output_filename(AMRScanFiles.Stage1.README_SUFFIX, config.card_version, ".txt"))
    try:
        with open(readme_path, 'w') as readme_file:
            readme_file.write("AMRScan_DB_CARD_prep Output Files Documentation\n")
            readme_file.write("==========================================\n\n")
            readme_file.write(f"CARD Version: {config.card_version}\n\n")
            readme_file.write("This document describes the final output files.\n\n")
            readme_file.write("--- Final FASTA Files ---\n")
            readme_file.write("Two final, non-redundant FASTA files are generated. They are perfectly synchronized; the Nth sequence in one file corresponds to the Nth sequence in the other, and they share the exact same FASTA header.\n\n")
            readme_file.write(f"1. **{format_output_filename(AMRScanFiles.Stage2.NR_FASTA, config.card_version, '.fasta')}**\n")
            readme_file.write("   - Contains the final, unique protein sequences and rRNA (nucleotide) sequences.\n\n")
            readme_file.write(f"2. **{format_output_filename(AMRScanFiles.Stage2.NR_FASTA_NUC, config.card_version, '.fasta')}**\n")
            readme_file.write("   - Contains the corresponding nucleotide sequences for every entry in the main FASTA file.\n\n")
            readme_file.write("**Header Format**\n")
            readme_file.write("Format: `>ARO_[AROs]__[GeneName]__[ModelCode][__Mutations]`\n")
            readme_file.write("Example: `>ARO_3000597__AAC_6_-Iad__V__S181G` (with mutation)\n")
            readme_file.write("Example: `>ARO_3000001__blaTEM-1__H` (wildtype, no mutation part)\n")
            readme_file.write(" - `[AROs]`: Consolidated ARO numbers for this sequence.\n")
            readme_file.write(" - `[GeneName]`: The CARD short name of a representative entry.\n")
            readme_file.write(" - `[ModelCode]`: Single character for model type (H:homolog, K:knockout, V:variant, O:overexpression, R:rRNA).\n")
            readme_file.write(" - `[__Mutations]`: Optional. Consolidated mutations, if present. Omitted for wildtype sequences.\n\n")
            readme_file.write(f"--- Final Metadata File: {format_output_filename(AMRScanFiles.Stage2.NR_METADATA, config.card_version, '.txt')} ---\n")
            readme_file.write("This file contains consolidated metadata for each unique sequence.\n\n")
            readme_file.write("Columns:\n\n")
            readme_file.write("1.  **Sequence_ID**: The unique ID that matches the FASTA headers. This is the primary key.\n")
            readme_file.write("2-8. **ARO_Number, etc.**: Metadata from a single representative pre-consolidation entry.\n")
            readme_file.write("9.  **Consolidated_Stage1_IDs**: A semicolon-separated list of all original IDs from Stage 1 that were merged into this single unique sequence. Provides a full audit trail.\n")
            readme_file.write("10. **Consolidated_ARO_Numbers**: All unique ARO numbers for this sequence.\n")
            readme_file.write("11. **Consolidated_Mutation_Strings**: All unique mutations that result in this sequence (empty if wildtype).\n")
            readme_file.write("12-17. Additional consolidated metadata (Parameter Types, Redundancy Count, Gene Family, etc.).\n")
        print(BColors.green(f"README file written to: {readme_path}"))
    except IOError as e:
        print(BColors.red(f"Error writing README file: {e}"), file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description="AMRScan_DB_CARD_prep - Processes CARD data into a non-redundant, annotated database."
    )
    parser.add_argument(
        "-c", "--config",
        required=True,
        help="Path to the configuration file (mandatory)."
    )
    parser.add_argument(
        "-d", "--output-dir",
        required=True,
        help="Name of the output directory (mandatory)."
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete and recreate output directory if it exists."
    )
    
    args = parser.parse_args()

    try:
        config = read_config(args.config)
        config.output_dir_name = args.output_dir
        config.overwrite_outputs = args.overwrite

        setup_output_directory(config.output_dir_name, config.overwrite_outputs)
        
        db_dir = config.output_dir_name
        # log_file_path = os.path.join(db_dir, format_output_filename(AMRScanFiles.Stage1.LOG_SUFFIX, config.card_version, ".log")) # Removed
        
        # with StreamRedirector(sys.stdout, sys.stderr, log_file_path): # Removed
        print(BColors.cyan(f"Starting AMRScan_DB_CARD_prep with CARD_VERSION: {config.card_version}"))
        
        sequence_to_nucleotide_map = {} # Stores protein_seq -> nucleotide_seq for Stage 2
        
        run_first_stage(config, db_dir, sequence_to_nucleotide_map)
        run_second_stage(config, db_dir, sequence_to_nucleotide_map)
        write_readme_file(config, db_dir)
        
        print(BColors.green("\n--- All Stages Complete ---"))
        print(BColors.green(f"Final processed files are in '{config.output_dir_name}' directory."))

    except Exception as e:
        sys.stderr.write(BColors.red(f"Critical Error: {e}\n"))
        sys.exit(1)

if __name__ == "__main__":
    main()
