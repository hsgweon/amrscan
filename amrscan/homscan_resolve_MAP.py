#!/usr/bin/env python3

# amrscan_pipeline/amrscan/homscan_resolve_MAP.py
import argparse
import sys
import re
import json
import csv
from collections import defaultdict

import pandas as pd
import numpy as np

pd.set_option('display.width', 120)
pd.set_option('display.max_columns', 15)
pd.set_option('display.max_colwidth', 40)

class BColors:
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    enabled = sys.stdout.isatty()
    @classmethod
    def _colorize(cls, c, t): return f"{c}{t}{cls.ENDC}" if cls.enabled else t
    @classmethod
    def cyan(cls, t): return cls._colorize(cls.OKCYAN, t)
    @classmethod
    def green(cls, t): return cls._colorize(cls.OKGREEN, t)
    @classmethod
    def red(cls, t): return cls._colorize(cls.FAIL, t)
    @classmethod
    def yellow(cls, t): return cls._colorize(cls.WARNING, t)

def parse_input_table(filepath: str, metric_col: str):
    print(BColors.cyan(f"--- Parsing input data from: {filepath} ---"))
    try:
        df = pd.read_csv(filepath, sep='\t')
    except FileNotFoundError:
        print(BColors.red(f"Error: Input file not found at '{filepath}'"), file=sys.stderr)
        sys.exit(1)
    if metric_col not in df.columns:
        print(BColors.red(f"Error: Metric column '{metric_col}' not found."), file=sys.stderr)
        sys.exit(1)
    
    df[metric_col] = pd.to_numeric(df[metric_col], errors='coerce').fillna(0)

    all_aros = set()
    aro_pattern = re.compile(r'ARO_\d+')
    family_map = {}
    
    for key in df['#AMR_Gene_Family;ARO']:
        parts = key.split(';', 1)
        family_part = parts[0]
        aros_in_key = aro_pattern.findall(key)
        all_aros.update(aros_in_key)
        for aro in aros_in_key:
            if aro not in family_map:
                family_map[aro] = family_part

    all_aros_list = sorted(list(all_aros))
    
    print(BColors.green(f"--- Found {len(df)} total observations for {len(all_aros_list)} unique AROs."))
    return df, all_aros_list, family_map

def parse_priors_file(filepath: str):
    if not filepath: return {}
    print(BColors.cyan(f"--- Loading priors from: {filepath} ---"))
    priors = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                parts = line.strip().split('\t')
                if len(parts) >= 1:
                    priors[parts[0]] = float(parts[1]) if len(parts) > 1 else 1.0
    except FileNotFoundError:
        print(BColors.yellow(f"Warning: Priors file not found at '{filepath}'. Proceeding without explicit priors."), file=sys.stderr)
        return {}
    print(BColors.green(f"--- Loaded {len(priors)} historical/clinical priors."))
    return priors

def run_iterative_solver(detailed_df, all_aros_list, priors, base_prior, prior_strength, metric_col: str, resolved_metric_col: str):
    """
    Resolves ambiguous abundances using a fast, iterative allocation algorithm (EM-like).
    This finds the Maximum A Posteriori (MAP) estimate.
    """
    print(BColors.cyan("--- Defining fast, iterative allocation solver... ---"))
    aro_pattern = re.compile(r'ARO_\d+')

    unique_abundances = defaultdict(float)
    for _, row in detailed_df.iterrows():
        key = row['#AMR_Gene_Family;ARO']
        if 'multiple[' not in key and not key.endswith(';multiple'):
            aros = aro_pattern.findall(key)
            if aros:
                unique_abundances[aros[0]] += row[metric_col]
    
    static_abundances = pd.Series(unique_abundances, dtype=float).reindex(all_aros_list, fill_value=0.0)

    for aro, hist_count in priors.items():
        if aro in static_abundances.index:
            static_abundances[aro] += hist_count * prior_strength * 0.1 
    
    static_abundances += base_prior * 1e-6

    ambiguous_rows = detailed_df[
        detailed_df['#AMR_Gene_Family;ARO'].str.contains('multiple', na=False)
    ].copy()
    
    ambiguous_rows['AROs_in_group'] = ambiguous_rows['#AMR_Gene_Family;ARO'].apply(
        lambda x: sorted(list(set(aro_pattern.findall(x))))
    )

    final_abundances = static_abundances.copy()

    print(BColors.cyan("--- Starting iterative allocation... ---"))
    for i in range(50):
        abundances_before_iter = final_abundances.copy()
        
        ambiguous_abun_this_iter = pd.Series(0.0, index=all_aros_list)

        for _, row in ambiguous_rows.iterrows():
            aros_in_group = row['AROs_in_group']
            if not aros_in_group: continue

            current_abun_of_group_members = abundances_before_iter.loc[aros_in_group]
            total_abun_for_group = current_abun_of_group_members.sum()

            if total_abun_for_group > 1e-9:
                proportions = current_abun_of_group_members / total_abun_for_group
            else:
                num_members = len(aros_in_group)
                if num_members > 0:
                    proportions = pd.Series(1.0 / num_members, index=aros_in_group)
                else:
                    continue

            distributed_abun = row[metric_col] * proportions
            ambiguous_abun_this_iter = ambiguous_abun_this_iter.add(distributed_abun, fill_value=0)

        final_abundances = static_abundances.add(ambiguous_abun_this_iter, fill_value=0)
        
        if np.allclose(final_abundances.values, abundances_before_iter.values, atol=1e-6, rtol=1e-5):
            print(BColors.green(f"--- Converged after {i+1} iterations. ---"))
            break
    else:
        print(BColors.yellow("--- Reached maximum iterations without convergence. ---"))

    result_df = final_abundances.reset_index()
    result_df.columns = ['ARO', resolved_metric_col]
    return result_df

def generate_final_summary(resolved_df, detailed_df, metric_col: str, resolved_metric_col: str):
    """
    Generates the final summary with user-friendly allocation proportions.
    """
    print(BColors.cyan("--- Generating final summary report... ---"))
    
    aro_pattern = re.compile(r'ARO_\d+')

    family_to_ambiguous_aros = defaultdict(set)
    for _, row in detailed_df.iterrows():
        if 'multiple[' in row['#AMR_Gene_Family;ARO'] or row['#AMR_Gene_Family;ARO'].endswith(';multiple'):
            family = row['#AMR_Gene_Family;ARO'].split(';')[0]
            aros = aro_pattern.findall(row['#AMR_Gene_Family;ARO'])
            family_to_ambiguous_aros[family].update(aros)

    summary_rows = []
    processed_families = set()

    for index, row in detailed_df.iterrows():
        key = row['#AMR_Gene_Family;ARO']
        family = key.split(';')[0]

        if 'multiple[' not in key and not key.endswith(';multiple'):
            new_row = row.to_dict()
            aros = aro_pattern.findall(key)
            if aros:
                aro = aros[0]
                new_row['Top_ARO'] = aro
                new_row['Allocation_Proportions'] = f"{aro}:100.0%"
                summary_rows.append(new_row)

        elif family not in processed_families:
            ambiguous_aros_for_family = list(family_to_ambiguous_aros[family])
            
            agg_row = {col: 'multiple' for col in detailed_df.columns}
            agg_row['#AMR_Gene_Family;ARO'] = f"{family};multiple"
            
            family_multiple_rows = detailed_df[
                (detailed_df['#AMR_Gene_Family;ARO'].str.startswith(f"{family};multiple"))
            ]
            
            cols_to_sum = ['Read_Count', 'Fragment_Count']
            all_numeric_cols = family_multiple_rows.select_dtypes(include=np.number).columns
            cols_to_average = [col for col in all_numeric_cols if col not in cols_to_sum]

            summed_vals = family_multiple_rows[cols_to_sum].sum().to_dict()
            averaged_vals = family_multiple_rows[cols_to_average].mean().to_dict()
            
            agg_row.update(summed_vals)
            agg_row.update(averaged_vals)

            ambiguous_results_df = resolved_df[resolved_df['ARO'].isin(ambiguous_aros_for_family)]
            total_resolved_abun_for_family = ambiguous_results_df[resolved_metric_col].sum()
            
            if not ambiguous_results_df.empty and total_resolved_abun_for_family > 1e-9:
                proportions = {
                    r['ARO']: r[resolved_metric_col] / total_resolved_abun_for_family
                    for _, r in ambiguous_results_df.iterrows() if r[resolved_metric_col] > 0
                }
                
                proportion_values = list(proportions.values())
                is_tie = len(proportion_values) > 1 and all(abs(p - proportion_values[0]) < 1e-5 for p in proportion_values)

                if is_tie:
                    agg_row['Top_ARO'] = 'N/A'
                else:
                    winner_row = ambiguous_results_df.loc[ambiguous_results_df[resolved_metric_col].idxmax()]
                    agg_row['Top_ARO'] = winner_row['ARO']
                
                sorted_proportions = sorted(proportions.items(), key=lambda item: item[1], reverse=True)
                
                proportions_str = ",".join([f"{aro}:{prop*100:.1f}%" for aro, prop in sorted_proportions])
                agg_row['Allocation_Proportions'] = proportions_str
            else:
                agg_row['Top_ARO'] = "None"
                agg_row['Allocation_Proportions'] = ""
            
            summary_rows.append(agg_row)
            processed_families.add(family)

    summary_df = pd.DataFrame(summary_rows)
    
    original_sum = detailed_df[metric_col].sum()
    resolved_sum = resolved_df[resolved_metric_col].sum()
    print(BColors.cyan("\n--- Conservation of Abundance Check (Internal Algorithm Validation) ---"))
    print(f"Sum of original '{metric_col}' column: {original_sum:.6f}")
    print(f"Sum of ALL individual resolved abundances (internal check): {resolved_sum:.6f}")
    if abs(original_sum - resolved_sum) < 1e-3:
        print(BColors.green("SUCCESS: Total abundance was conserved by the solver."))
    else:
        print(BColors.red(f"ERROR: Total abundance was NOT conserved by the solver. Difference: {original_sum - resolved_sum:.6f}"))

    return summary_df

def main():
    parser = argparse.ArgumentParser(
        description="From homscan data, resolves ambiguous gene families to determine the most likely top gene (Top_ARO) and the proportion of evidence allocated to each candidate gene. Does not create a new abundance column, but reports proportions.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-i", "--input-file", required=True, help="Path to the input `_homscan_detailed.tsv` file.")
    parser.add_argument("-p", "--priors-file", help="Optional: Path to a tab-separated file of priors.")
    parser.add_argument("-o", "--output-prefix", required=True, help="Prefix for the output file (e.g., 'SRR12345').")
    parser.add_argument("--metric-column", default="RPKG", help="The numeric column to use as the abundance metric (default: RPKG).")
    parser.add_argument("--base-prior", type=float, default=1.0, help="Baseline prior 'pseudo-count' for genes NOT in the priors file (default: 1.0).")
    parser.add_argument("--prior-strength", type=float, default=1.0, help="Multiplier for the influence of the priors file (default: 1.0).")
    args = parser.parse_args()

    output_filename = f"{args.output_prefix}_homscan_MAP.tsv"
    resolved_metric_col = f"Resolved_{args.metric_column}"

    detailed_df, all_aros_list, family_map = parse_input_table(args.input_file, args.metric_column)
    priors = parse_priors_file(args.priors_file)
    
    if detailed_df.empty or detailed_df[args.metric_column].sum() == 0:
        print(BColors.yellow("Warning: No valid data or abundance found in the input file. Writing an empty output."), file=sys.stderr)
        with open(output_filename, 'w') as f:
            f.write("")
        sys.exit(0)

    resolved_df = run_iterative_solver(detailed_df, all_aros_list, priors, args.base_prior, args.prior_strength, args.metric_column, resolved_metric_col)
    summary_df = generate_final_summary(resolved_df, detailed_df, args.metric_column, resolved_metric_col)
    
    try:
        homscan_cols = [
            '#AMR_Gene_Family;ARO', 'Read_Count', 'Fragment_Count', 
            'Lateral_Coverage_%', 'Gene_Length_bp', 'RPK', 'FPK', 'RPKG', 'FPKG', 
            'RPKPC', 'FPKPC', 'RPKPMC', 'FPKPMC'
        ]
        
        for col in homscan_cols:
            if col not in summary_df.columns:
                summary_df[col] = 0

        new_cols = ['Top_ARO', 'Allocation_Proportions']
        final_cols = homscan_cols + new_cols
        
        summary_df = summary_df.reindex(columns=final_cols)
        
        summary_df.to_csv(output_filename, sep='\t', index=False, float_format='%.6f', quoting=csv.QUOTE_NONE, escapechar='\\')
        print(BColors.green(f"\n--- Successfully wrote final summary to: {output_filename} ---"))
        
        print("\n--- Preview of Final Summary ---")
        preview_cols = ['#AMR_Gene_Family;ARO', args.metric_column, 'Top_ARO', 'Allocation_Proportions']
        pd.set_option('display.max_colwidth', 80)
        print(summary_df[preview_cols].head(20).to_string(index=False))

    except IOError as e:
        print(BColors.red(f"Error writing output file: {e}"), file=sys.stderr)
        sys.exit(1)
    except KeyError as e:
        print(BColors.red(f"Error: A required column was not found during final processing: {e}"), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()