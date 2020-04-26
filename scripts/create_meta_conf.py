#!/usr/bin/env python3

import argparse
import json

def run():

    parser = argparse.ArgumentParser(description='Create a json config file for meta-analysis')
    parser.add_argument('in_files_loc', action='store', type=str, help='A tab-delimited text file with gs:// locations of sumstat input files. One row is one meta')
    parser.add_argument('out_file_prefix', action='store', type=str, help='Output json file prefix')
    args = parser.parse_args()

    with open(args.in_files_loc) as f:
        index = 1
        for line in f:
            conf = {'meta': []}
            for file in line.strip().split('\t'):
                # files = [file for file in [line.strip().split('\t') for line in lines]]
                # files = [file for sublist in files for file in sublist]
                filename = file.split('/').pop()
                fields = filename.replace('.formatted', '').replace('.munged', '').replace('.gz', '').replace('.txt', '').split('.')
                if (len(fields) != 10):
                    raise Exception('Unexpected filename in {}: {}'.format(args.in_files_loc, filename))
                try:
                    n_cases = int(fields[6])
                    n_controls = int(fields[7])
                except ValueError:
                    raise Exception('Could not parse numbers of cases and controls from {}'.format(filename))
                conf['meta'].append({
                    'name': fields[0],
                    'file': file.replace('gs://', '/cromwell_root/'),
                    'n_cases': n_cases,
                    'n_controls': n_controls,
        	        'chr':'#CHR',
                    'pos':'POS',
                    'ref':'Allele1',
                    'alt':'Allele2',
                    'effect':'BETA',
                    'effect_type':'beta',
                    'pval':'p.value',
                    'se':'SE',
                    'extra_cols':['AF_Allele2', 'imputationInfo']
                })
            outfile = args.out_file_prefix + '_' + str(index) + '.json'
            with open(outfile, 'w') as out:
                json.dump(conf, out, indent=4)
                print(outfile + ' written')
            index = index + 1

if __name__ == '__main__':
    run()
