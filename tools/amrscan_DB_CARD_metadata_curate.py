#!/usr/bin/env python3

import argparse
import os
import pandas as pd

class BColors:
    """A helper class to add color to terminal output."""
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    enabled = True

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

def load_rules(rules_path):
    """Loads curation rules from a tab-separated file."""
    print(BColors.cyan(f"--- Loading curation rules from: {rules_path} ---"))
    if not os.path.exists(rules_path):
        print(BColors.red(f"Error: Rules file not found at '{rules_path}'"), file=sys.stderr)
        return None
    
    rules = []
    try:
        with open(rules_path, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        rules.append((parts[0].strip(), parts[1].strip()))
                    else:
                        print(BColors.yellow(f"Warning: Skipping malformed rule line: {line.strip()}"))
    except IOError as e:
        print(BColors.red(f"Error reading rules file: {e}"), file=sys.stderr)
        return None
    
    print(BColors.green(f"--- Loaded {len(rules)} curation rules."))
    return rules

def main():
    parser = argparse.ArgumentParser(
        description="Curates an AMRScan metadata file by applying custom renaming rules to the AMR_Gene_Family column."
    )
    parser.add_argument(
        "-m", "--metadata",
        required=True,
        help="Path to the AMRScan metadata file to be curated (e.g., amrscan_DB_CARD_v4.0.0_NR_metadata.txt)."
    )
    parser.add_argument(
        "-r", "--rules",
        required=True,
        help="Path to a tab-separated text file with curation rules. Format: 'ARO_Name_Contains\\tNew_AMR_Gene_Family'."
    )
    parser.add_argument(
        "-o", "--output",
        help="Path for the new, curated metadata file. If not provided, the original file will be overwritten (in-place)."
    )
    args = parser.parse_args()

    # 1. Load the rules
    rules = load_rules(args.rules)
    if rules is None:
        sys.exit(1)

    # 2. Load the metadata file into a pandas DataFrame
    print(BColors.cyan(f"--- Loading metadata from: {args.metadata} ---"))
    try:
        df = pd.read_csv(args.metadata, sep='\t')
    except FileNotFoundError:
        print(BColors.red(f"Error: Metadata file not found at '{args.metadata}'"), file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(BColors.red(f"Error parsing metadata file: {e}"), file=sys.stderr)
        sys.exit(1)

    # 3. Apply the curation rules
    print(BColors.cyan("--- Applying curation rules... ---"))
    changes_made = 0
    for index, row in df.iterrows():
        sequence_id = row['Sequence_ID']
        try:
            # Extract ARO_Name from format like 'ARO_3001770__OXA-43__H'
            aro_name = sequence_id.split('__')[1]
        except IndexError:
            continue # Skip if the header format is unexpected

        original_family = row['AMR_Gene_Family']
        new_family = original_family

        # Iterate through rules and apply the first one that matches
        for name_contains, replacement_family in rules:
            if name_contains in aro_name:
                new_family = replacement_family
                break # Stop after the first matching rule
        
        if new_family != original_family:
            df.at[index, 'AMR_Gene_Family'] = new_family
            changes_made += 1
            print(f"Changed '{sequence_id}': '{original_family}' -> '{new_family}'")

    print(BColors.green(f"--- Curation complete. Made {changes_made} changes. ---"))

    # 4. Save the result
    output_path = args.output if args.output else args.metadata
    print(BColors.cyan(f"--- Saving curated metadata to: {output_path} ---"))
    try:
        df.to_csv(output_path, sep='\t', index=False)
        print(BColors.green("--- Successfully saved file."))
    except IOError as e:
        print(BColors.red(f"Error writing output file: {e}"), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()