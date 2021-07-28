#!/usr/bin/env python3

import argparse
import json
import pandas as pd


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
        'extra_cols': ['af_alt', 'af_alt_cases', 'af_alt_controls', 'nearest_genes']
    }

    ukbb_conf = {
        'name': 'UKBB',
        'file': mapping['ukbb_link'],
        'n_cases': int(float(mapping['ukbb_n_cases'])),
        'n_controls': int(float(mapping['ukbb_n_controls'])),
        'chr': '#chr',
        'pos': 'pos',
        'ref': 'ref',
        'alt': 'alt',
        'effect': 'beta_EUR',
        'effect_type': 'beta',
        'pval': 'pval_EUR',
        'se': 'se_EUR',
        'extra_cols': ['low_confidence_EUR']
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
        'extra_cols': ["af_alt", "AF.Cases", "AF.Controls"]
    }

    with open(mapping['fg_phenotype'] + '.json', 'w') as out:
        json.dump({'meta': [fg_conf, ukbb_conf, est_conf]}, out, indent=4)


def run():

    parser = argparse.ArgumentParser(description='Create a json config file for meta-analysis combining results from FinnGen, UKBB and EST')
    parser.add_argument('in_mapping_file', action='store', type=str, help='A tab-delimited text file with phenotype mapping information between studies')

    args = parser.parse_args()

    # Read mapping and drop rows with missing values (missing pair)
    mapping = pd.read_csv(args.in_mapping_file, sep='\t', dtype=object)
    mapping.dropna(axis=0, how='any', inplace=True)
    #mapping.index = pd.MultiIndex.from_frame(mapping[['fg_phenotype', 'ukbb_phenotype', 'est_phenotype']])

    # Print summary stat links to file
    mapping[['fg_link', 'ukbb_link', 'estbb_link']].to_csv('sumstat_files.txt', header=False, index=False, sep='\t')

    # Replace gs prefix with cromwell_root for cromwell run
    mapping['fg_link'] = mapping['fg_link'].replace('^gs:/', '/cromwell_root', regex=True)
    mapping['ukbb_link'] = mapping['ukbb_link'].replace('^gs:/', '/cromwell_root', regex=True)
    mapping['estbb_link'] = mapping['estbb_link'].replace('^gs:/', '/cromwell_root', regex=True)

    # Generate json configs for meta-analysis
    mapping.apply(generate_json, axis=1)


if __name__ == '__main__':
    run()
