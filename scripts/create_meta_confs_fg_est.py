#!/usr/bin/env python3

import argparse
import json
import pandas as pd
import os
import subprocess
import shlex

REQUIRED_COLS = ['fg_phenotype', 'fg_link', 'estbb_link']
OPTIONAL_COLS = ['fg_n_cases', 'estbb_n_cases', 'fg_n_controls', 'estbb_n_controls']


def generate_json(mapping):

    fg_conf = {
        'name': 'FINNGEN',
        'file': mapping['fg_link'],
        'n_cases': int(float(mapping['fg_n_cases'])),
        'n_controls': int(float(mapping['fg_n_controls'])),
        'chr': '#chrom',
        'pos': 'pos',
        'ref': 'ref',
        'alt': 'alt',
        'effect': 'beta',
        'effect_type': 'beta',
        'pval': 'pval',
        'se': 'sebeta',
        'extra_cols': ['af_alt', 'af_alt_cases', 'af_alt_controls']
    }

    est_conf = {
        'name': 'ESTBB',
        'file': mapping['estbb_link'],
        'n_cases': int(float(mapping['estbb_n_cases'])),
        'n_controls': int(float(mapping['estbb_n_controls'])),
        'chr': '#CHR',
        'pos': 'POS',
        'ref': 'REF',
        'alt': 'ALT',
        'effect': 'beta',
        'effect_type': 'beta',
        'pval': 'pval',
        'se': 'sebeta',
        'extra_cols': ["af_alt"]
    }

    json_file = 'jsons/' + mapping['fg_phenotype'] + '.json' if os.path.isdir('jsons') else mapping['fg_phenotype'] + '.json'
    with open(json_file, 'w') as out:
        json.dump({'meta': [fg_conf, est_conf]}, out, indent=4)


def run():

    parser = argparse.ArgumentParser(description='Create a json config file for meta-analysis combining results from FinnGen and EST')
    parser.add_argument('in_mapping_file', action='store', type=str, help='A tab-delimited text file with phenotype mapping information between studies')
    parser.add_argument('--bucket', action='store', type=str, help='GCS bucket path')
    parser.add_argument('--sumstat_filelist_name', action='store', type=str, default='sumstat_files.txt', help='Name of file listing the sumstats')
    parser.add_argument('--json_filelist_name', action='store', type=str, default='conf_jsons.txt', help='Name of file listing the jsons')

    args = parser.parse_args()

    # Check mapping is ok and filter
    print('Checking mapping...')
    mapping = pd.read_csv(args.in_mapping_file, sep='\t', dtype=object)
    missing_cols = [col for col in REQUIRED_COLS if col not in mapping.columns]
    if len(missing_cols) > 0:
        raise Exception('Missing required columns: ' + ', '.join(missing_cols))
    for col in OPTIONAL_COLS:
        if col not in mapping.columns:
            mapping[col] = 0
    mapping = mapping[REQUIRED_COLS + OPTIONAL_COLS]
    mapping[OPTIONAL_COLS] = mapping[OPTIONAL_COLS].fillna(0)
    mapping.dropna(axis=0, how='any', inplace=True)

    # Check for duplicate phenotype names
    if any(mapping[['fg_phenotype']].duplicated()):
        duplicates = pd.unique(mapping.loc[mapping[['fg_phenotype']].duplicated(), 'fg_phenotype']).tolist()
        raise Exception(f'Found duplicate phenotype names: {duplicates}.\n Make them unique and try again...')

    # Print summary stat links to file
    mapping[['fg_link', 'estbb_link']].to_csv(args.sumstat_filelist_name, header=False, index=False, sep='\t')

    # Replace gs prefix with cromwell_root for cromwell run
    mapping['fg_link'] = mapping['fg_link'].replace('^gs:/', '/cromwell_root', regex=True)
    mapping['estbb_link'] = mapping['estbb_link'].replace('^gs:/', '/cromwell_root', regex=True)

    # Generate json configs for meta-analysis
    print('Generating config files...')
    try:
        os.mkdir('jsons')
    except OSError as error:
        print(error)
    mapping.apply(generate_json, axis=1)

    if args.bucket is not None:
        with open(args.json_filelist_name, 'w') as f:
            for pheno in mapping['fg_phenotype']:
                f.write(os.path.join(args.bucket, 'jsons', pheno) + '.json\n')
        print(f'Transferring files to bucket {args.bucket} ...')
        cmd_1 = f'gsutil -q cp {args.sumstat_filelist_name} {args.bucket}'
        cmd_2 = f'gsutil -q cp {args.json_filelist_name} {args.bucket}'
        cmd_3 = f'gsutil -q -m cp -r jsons/ {args.bucket}'
        subprocess.run(shlex.split(cmd_1))
        subprocess.run(shlex.split(cmd_2))
        subprocess.run(shlex.split(cmd_3))
    
    print('Done.')


if __name__ == '__main__':
    run()
