#!/usr/bin/env python

import sys
import math
import argparse
import gzip

# Example usage:
# python3 munge.py biomarkers-30780-both_sexes-irnt.tsv.bgz --chr-col chr --pos-col pos --ref-col ref --alt-col alt --af-col af_EUR --effect-col beta_EUR --se-col se_EUR --pval-col neglog10_pval_EUR --effect-type beta --pval-type mlog10p | sort -k 1,1g -k2,2g | bgzip > munged.gz

def process_file(f_in, args):
    """
    Reads and processes the summary statistics file line by line.
    """
    # Normalize settings from args
    eff_type = args.effect_type.lower()
    s_type = args.se_type.lower()
    p_type = args.pval_type.lower()
    flip = args.flip_alleles
    run_filter = True if args.filt_col and args.filt_threshold is not None else False
    af_allele_is_ref = (args.af_allele.lower() == "ref")
    filt_threshold = args.filt_threshold

    # Create the dynamic mapping of input column names to our standard names
    col_rename_map = {
        args.chr_col: "#CHR",
        args.pos_col: "POS",
        args.ref_col: "REF",
        args.alt_col: "ALT",
        args.af_col: "af_alt",
        args.effect_col: "beta",
        args.se_col: "sebeta",
        args.pval_col: "pval",
    }

    col_indices = {}
    pos_idx = -1
    extra_idx = -1
    filt_idx = -1
    col_map_name_to_index = {} # Will map {standard_name: index}
    error_lines = 0

    # Process every line from the input file handle
    for line_num, line in enumerate(f_in):
        try:
            # Split line on any whitespace (tabs or spaces)
            fields = line.strip().split(args.delim)
            if not fields:
                continue  # Skip empty lines

            # --- HEADER PROCESSING ---
            if line_num == 0:
                processed_header = []

                for i, col_name in enumerate(fields):
                    # Rename the column if it's in our map, otherwise keep the original name
                    standard_name = col_rename_map.get(col_name, col_name)
                    processed_header.append(standard_name)
                    # Store mapping for *all* columns (standardized or not)
                    col_map_name_to_index[standard_name] = i
                
                # Add the mlogp column to the header
                processed_header.append("mlogp")

                # Print the new header, joined by tabs (OFS="\t")
                print("\t".join(processed_header))

                # Now, get the 0-based indices for all columns we need.
                # Wrap this in a try/except for user-friendly error messages.
                try:
                    col_indices = {
                        "chr": col_map_name_to_index["#CHR"],
                        "pos": col_map_name_to_index["POS"],
                        "ref": col_map_name_to_index["REF"],
                        "alt": col_map_name_to_index["ALT"],
                        "af": col_map_name_to_index["af_alt"],
                        "beta": col_map_name_to_index["beta"],
                        "se": col_map_name_to_index["sebeta"],
                        "pval": col_map_name_to_index["pval"],
                    }
                    pos_idx = col_indices["pos"]
                    extra_idx = col_map_name_to_index.get("EXTRA", -1) 
                    
                    if run_filter:
                        # Check if the user-specified filter column exists
                        filt_idx = col_map_name_to_index[args.filt_col]
                        
                except KeyError as e:
                    sys.stderr.write(f"\nError: Missing required column.\n")
                    sys.stderr.write(f"Could not find standardized column {e} in the header.\n")
                    sys.stderr.write("Please check that your column name parameters (e.g., --chr-col, --pval-col) match the input file header.\n\n")
                    sys.exit(1)
                
                continue  # Done processing header, move to next line

            # --- DATA LINE PROCESSING ---

            # Filtering
            if run_filter:
                filt_val = float(fields[filt_idx])
                if filt_val <= filt_threshold:
                    continue  # Filter failed, skip this line

            # Get key values that need transformation
            beta_str = fields[col_indices["beta"]]
            se_str = fields[col_indices["se"]]
            pval_str = fields[col_indices["pval"]]
            af_str = fields[col_indices["af"]]

            # Effect type transformation (OR -> Beta)
            beta = float(beta_str)
            if eff_type == "or":
                beta = math.log(beta)

            # Standard Error transformation (CI -> SE)
            if s_type == "ci":
                ci_lower_str, ci_upper_str = se_str.split(',')
                ci_lower = float(ci_lower_str)
                ci_upper = float(ci_upper_str)
                
                if eff_type == "or":
                    ci_lower = math.log(ci_lower)
                    ci_upper = math.log(ci_upper)
                
                se = (ci_upper - ci_lower) / (2 * 1.96)
            else:
                se = float(se_str) # Already an SE

            # P-value transformation (mlog10p -> P)
            pval_raw = float(pval_str)
            if p_type == "mlog10p":
                pval = 10 ** (-pval_raw)
                mlogp = pval_raw
            else:
                pval = pval_raw
                if pval > 0:
                    mlogp = -math.log10(pval)
                else:
                    mlogp = float('inf')
                
            af = float(af_str)

            # Chromosome formatting
            chr_idx = col_indices["chr"]
            chrom = fields[chr_idx]
            chrom = chrom.lstrip('0')
            chrom = chrom.lstrip('chr')
            if chrom == "X":
                chrom = "23"
            elif chrom == "Y":
                chrom = "24"
            fields[chr_idx] = chrom

            # Check for "TEST_FAIL"
            is_fail = False
            if extra_idx != -1 and fields[extra_idx] == "TEST_FAIL":
                is_fail = True

            # Main QC Filter check
            qc_pass = (
                chrom.isdigit() and
                0 <= pval < 1 and
                0 < mlogp and
                -1e6 < beta < 1e6 and beta != 0 and
                0 < af < 1 and
                se > 0 and
                not is_fail
            )

            if qc_pass:
                # Allele flip logic
                if flip:
                    ref_idx, alt_idx = col_indices["ref"], col_indices["alt"]
                    fields[ref_idx], fields[alt_idx] = fields[alt_idx], fields[ref_idx]
                    beta = -beta
                    af = 1.0 - af

                # Allele frequency adjustment
                if af_allele_is_ref:
                    af = 1.0 - af

                # Update the fields list with transformed values
                fields[col_indices["beta"]] = str(beta)
                fields[col_indices["se"]] = str(se)
                fields[col_indices["pval"]] = str(pval)
                fields[col_indices["af"]] = str(af)
                
                # Format POS as an integer
                fields[pos_idx] = str(int(float(fields[pos_idx])))

                # Append the mlogp value
                fields.append(str(mlogp))

                # Print the processed line
                print("\t".join(fields))

        except (ValueError, IndexError, TypeError, ZeroDivisionError) as e:
            # Skip any lines that cause conversion or processing errors
            error_lines += 1
            if args.verbose:
                sys.stderr.write(f"Skipping bad line {line_num+1}: {line.strip()}. Error: {e}\n")
            continue
    sys.stderr.write(f"Finished processing. Skipped {error_lines} lines due to errors.\n")
    

def main():
    parser = argparse.ArgumentParser(
        description="Munge a summary statistics file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # --- Positional Argument ---
    parser.add_argument(
        "sumstat_file", 
        help="Input summary statistics file. Can be plain text or gzipped (.gz)."
    )
    
    # --- Column Name Arguments ---
    col_group = parser.add_argument_group('Column Mapping (Required)')
    col_group.add_argument("--chr-col", required=True, help="Name of the chromosome column.")
    col_group.add_argument("--pos-col", required=True, help="Name of the position column.")
    col_group.add_argument("--ref-col", required=True, help="Name of the reference/non-effect allele column.")
    col_group.add_argument("--alt-col", required=True, help="Name of the alternate/effect allele column.")
    col_group.add_argument("--af-col", required=True, help="Name of the allele frequency column.")
    col_group.add_argument("--effect-col", required=True, help="Name of the effect size column (beta or OR).")
    col_group.add_argument("--se-col", required=True, help="Name of the standard error or confidence interval column.")
    col_group.add_argument("--pval-col", required=True, help="Name of the p-value column.")

    # --- Transformation and Filtering Settings ---
    set_group = parser.add_argument_group('Transformation and Filtering')
    set_group.add_argument("--effect-type", choices=['beta', 'or'], default='beta', help="Type of effect size reported.")
    set_group.add_argument("--se-type", choices=['se', 'ci'], default='se', help="Type of error reported (standard error or confidence interval).")
    set_group.add_argument("--pval-type", choices=['p', 'mlog10p'], default='p', help="Type of p-value reported (raw p-value or -log10(p)).")
    set_group.add_argument("--af-allele", choices=['alt', 'ref'], default='alt', help="Allele described by the frequency in --af-col.")
    set_group.add_argument("--flip-alleles", action='store_true', help="Set this flag to flip REF/ALT alleles and invert beta/AF.")
    
    set_group.add_argument("--filt-col", help="Name of the column to use for filtering.")
    set_group.add_argument("--filt-threshold", type=float, help="Rows with a value in --filt-col LESS than or EQUAL to this are removed.")

    # --- Other Settings ---
    parser.add_argument("--delim", default="\t", help="Delimiter used in the summary statistics file.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output for debugging.")

    args = parser.parse_args()

    # --- Input Validation ---
    if args.filt_col and not args.filt_threshold or args.filt_threshold and not args.filt_col:
        parser.error("Both --filt-col and --filt-threshold must be provided together to enable filtering.")

    # --- File Handling ---
    try:
        if args.sumstat_file.endswith('.gz') or args.sumstat_file.endswith('.bgz'):
            with gzip.open(args.sumstat_file, 'rt', encoding='utf-8') as f_in:
                process_file(f_in, args)
        else:
            with open(args.sumstat_file, 'r', encoding='utf-8') as f_in:
                process_file(f_in, args)
                
    except FileNotFoundError:
        sys.stderr.write(f"Error: Input file not found: {args.sumstat_file}\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"An unexpected error occurred: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()