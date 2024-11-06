#!/usr/bin/env python3

import argparse
import json
import os
import os.path
import subprocess
import shlex
import gzip

def format_int(val):
    try:
        return int(val)
    except ValueError:
        return 0

def generate_json(data, h_idx, columns, studies):

    confs = []

    for study in studies:
        study_conf = {
            "name": study,
            "file": data[h_idx[columns["studies"][study]["link"]]].replace("gs:/", "/cromwell_root"),
            "n_cases": format_int(data[h_idx[columns["studies"][study]["n_cases"]]]),
            "n_controls": format_int(data[h_idx[columns["studies"][study]["n_controls"]]]),
            "chr": columns["studies"][study]["chr"],
            "pos": columns["studies"][study]["pos"],
            "ref": columns["studies"][study]["ref"],
            "alt": columns["studies"][study]["alt"],
            "effect": columns["studies"][study]["effect"],
            "effect_type": columns["studies"][study]["effect_type"],
            "pval": columns["studies"][study]["pval"],
            "extra_cols": columns["studies"][study]["extra_cols"]
        }
        if "se" in columns["studies"][study]:
            study_conf["se"] = columns["studies"][study]["se"]
        confs.append(study_conf)

    json_file = "jsons/" + data[h_idx[columns["phenotype"]]] + ".json" if os.path.isdir("jsons") else data[h_idx[columns["phenotype"]]] + ".json"
    with open(json_file, "w") as out:
        json.dump({"meta": confs}, out, indent=4)

def get_studies(header, link_suffix):
    return [col[:-len(link_suffix)] for col in header if col.endswith(link_suffix)]

def categorize_columns(header, args):
    columns = {}
    if args.study_cols_json:
        with open(args.study_cols_json, "r") as f:
            columns["studies"] = json.load(f)
    else:
        for s in args.studies:
            columns["studies"][s] = {
                "chr": "#CHR",
                "pos": "POS",
                "ref": "REF",
                "alt": "ALT",
                "effect": "beta",
                "effect_type": "beta",
                "pval": "pval",
                "se": "sebeta",
                "extra_cols": ["af_alt"]
            }
    
    for col in header:
        if col == args.phenotype_col:
            columns["phenotype"] = col
        for s in args.studies:
            if col.startswith(s):
                if col.endswith(args.link_col_suffix):
                    columns["studies"][s]["link"] = col
                elif col.endswith(args.n_cases_col_suffix):
                    columns["studies"][s]["n_cases"] = col
                elif col.endswith(args.n_controls_col_suffix):
                    columns["studies"][s]["n_controls"] = col

    return columns


def run():

    parser = argparse.ArgumentParser(description="Create a json config file for meta-analysis")
    parser.add_argument("in_mapping_file", action="store", type=str, help="A tab-delimited text file with phenotype mapping information between studies")
    parser.add_argument("--studies", action="store", type=str, nargs="+", help="List of studies to include in the meta-analysis configs. (Default: all studies in the mapping file)")
    parser.add_argument("--bucket", action="store", type=str, help="GCS bucket path")
    parser.add_argument("--sumstat_filelist_name", action="store", type=str, default="sumstat_files.txt", help="Name of file listing the sumstats")
    parser.add_argument("--json_filelist_name", action="store", type=str, default="conf_jsons.txt", help="Name of file listing the jsons")
    parser.add_argument("--phenotype_col", action="store", type=str, default="phenotype", help="Column name containing the phenotype name. (Default: 'phenotype')")
    parser.add_argument("--link_col_suffix", action="store", type=str, default="_link", help="Suffix for the column name containing the sumstat URI. (Default: '_link')")
    parser.add_argument("--n_cases_col_suffix", action="store", type=str, default="_n_cases", help="Suffix for the column name containing the number of cases. (Default: '_n_cases')")
    parser.add_argument("--n_controls_col_suffix", action="store", type=str, default="_n_controls", help="Suffix for the column name containing the number of controls. (Default: '_n_controls')")
    parser.add_argument("--study_cols_json", action="store", type=str, help="JSON file containing the column names in sumstat files per study")
    parser.add_argument("--min_studies", action="store", type=int, default=2, help="Minimum number of studies required for a phenotype to be included in the meta-analysis. (Default: 2)")
    parser.add_argument("--complete", action="store_true", help="Only include phenotypes with all studies present")
    parser.add_argument("--required_studies", action="store", type=str, nargs="+", help="List of studies that are required for a phenotype to be included in the meta-analysis")

    args = parser.parse_args()

    # Check mapping is ok and filter
    print("Generating config files...")
    print("Skipping phenotypes with names as 'NA'...")

    try:
        os.mkdir("jsons")
    except OSError as error:
        print(error)

    with gzip.open(args.in_mapping_file, "rt") if args.in_mapping_file.endswith(".gz") else open(args.in_mapping_file, "r") as f_in, open(args.sumstat_filelist_name, "w") as f_out, open(args.json_filelist_name, "w") as f_out_json:
        header = f_in.readline().strip().split("\t")
        if not args.studies:
            args.studies = get_studies(header, args.link_col_suffix)
        if args.complete:
            args.min_studies = len(args.studies)
            print(f"Only including phenotypes with all studies ({args.min_studies}) present.")
        columns = categorize_columns(header, args)
        h_idx = {h:i for i,h in enumerate(header)}
        for line in f_in:
            d = line.strip().split("\t")
            if d[h_idx[columns["phenotype"]]] == "NA":
                continue
            links = [d[h_idx[columns["studies"][s]["link"]]] for s in columns["studies"] if d[h_idx[columns["studies"][s]["link"]]] != "NA"]
            pheno_studies = [s for s in columns["studies"] if d[h_idx[columns["studies"][s]["link"]]] != "NA"]
            if len(links) < args.min_studies:
                print(f"Skipping {d[h_idx[columns['phenotype']]]} as it has less than {args.min_studies} studies.")
                continue
            if args.required_studies:
                if not all([study in pheno_studies for study in args.required_studies]):
                    print(f"Skipping {d[h_idx[columns['phenotype']]]} as it does not have all required studies.")
                    continue
            f_out.write("\t".join(links) + "\n")
            generate_json(d, h_idx, columns, pheno_studies)
            if args.bucket is not None:
                f_out_json.write(os.path.join(args.bucket, "jsons", d[h_idx[columns["phenotype"]]]) + ".json\n")

    if args.bucket is not None:
        print(f"Transferring files to bucket {args.bucket} ...")
        cmd_1 = f"gsutil -q cp {args.sumstat_filelist_name} {args.bucket}"
        cmd_2 = f"gsutil -q cp {args.json_filelist_name} {args.bucket}"
        cmd_3 = f"gsutil -q -m cp -r jsons/ {args.bucket}"
        subprocess.run(shlex.split(cmd_1))
        subprocess.run(shlex.split(cmd_2))
        subprocess.run(shlex.split(cmd_3))
    
    print("Done.")


if __name__ == "__main__":
    run()
