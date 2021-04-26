#!/usr/bin/env python3

import argparse
import json
import re
import pandas as pd


def generate_json(combined_results):

    fg_conf = {
        'name': 'FINNGEN',
        'file': combined_results['fg_link'],
        'n_cases': int(combined_results['num_cases']),
        'n_controls': int(combined_results['num_controls']),
        'chr': '#chrom',
        'pos': 'pos',
        'ref': 'ref',
        'alt': 'alt',
        'effect': 'beta',
        'effect_type': 'beta',
        'pval': 'pval',
        'se': 'sebeta',
        'extra_cols': ['af_alt', 'af_alt_cases', 'af_alt_controls', 'rsids']
    }

    ukbb_conf = {
        'name': 'UKBB',
        'file': combined_results['ukbb_link'],
        'n_cases': int(combined_results['n_cases_EUR']),
        'n_controls': int(combined_results['n_controls_EUR']),
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

    with open(combined_results['fg_phenotype'] + '.json', 'w') as out:
        json.dump({'meta': [fg_conf, ukbb_conf]}, out, indent=4)


def run():

    parser = argparse.ArgumentParser(description='Create a json config file for meta-analysis combining results from FinnGen and UKBB')
    parser.add_argument('in_mapping_file', action='store', type=str, help='A tab-delimited text file with phenotype mapping information between studies')
    parser.add_argument('in_study_fg_file', action='store', type=str, help='A tab-delimited text file with FinnGen phenotype N\'s')
    parser.add_argument('in_study_ukbb_file', action='store', type=str, help='A tab-delimited text file with UKBB phenotype N\'s')
    args = parser.parse_args()

    # Read mapping and drop rows with missing values (missing pair)
    mapping = pd.read_csv(args.in_mapping_file, sep='\t', dtype=object)
    mapping.dropna(axis=0, how='any', inplace=True)
    mapping.index = pd.MultiIndex.from_frame(mapping[['fg_phenotype', 'ukbb_phenotype']])

    # Read fg study metadata and drop phenotypes not in mapping
    fg = pd.read_csv(args.in_study_fg_file, sep='\t')
    fg = fg[fg['phenocode'].isin(mapping['fg_phenotype'])]
    fg.index = pd.Index(fg.phenocode, name='fg_phenotype')

    # Read ukbb study metadata and drop phenotypes not in mapping
    # Note: UKBB has duplicate phenocodes.
    # Remove from duplicate pairs ones with 'trait_type' == 'continuous'
    ukbb = pd.read_csv(args.in_study_ukbb_file, sep='\t')
    ukbb = ukbb[ukbb['phenocode'].isin(mapping['ukbb_phenotype'])]
    duplicated_ids = ukbb['phenocode'][ukbb['phenocode'].duplicated()]
    ukbb = ukbb[~((ukbb['phenocode'].isin(duplicated_ids)) & (ukbb['trait_type'] == 'continuous'))]
    ukbb.index = pd.Index(ukbb.phenocode, name='ukbb_phenotype')

    # Join study metadatas with mapping
    all = mapping.join(fg, how='inner').join(ukbb, how='inner', rsuffix='_ukbb')

    # Print summary stat links to file
    #all['fg_link'].to_csv('fg_files.txt', header=False, index=False)
    #all['ukbb_link'].to_csv('ukbb_files.txt', header=False, index=False)
    all[['fg_link', 'ukbb_link']].to_csv('sumstat_files.txt', header=False, index=False, sep='\t')

    # Replace gs prefix with cromwell_root for cromwell run
    all['fg_link'] = all['fg_link'].replace('^gs:/', '/cromwell_root', regex=True)
    all['ukbb_link'] = all['ukbb_link'].replace('^gs:/', '/cromwell_root', regex=True)

    # Generate json configs for meta-analysis
    all.apply(generate_json, axis=1)


if __name__ == '__main__':
    run()
