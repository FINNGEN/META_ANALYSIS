#!/usr/bin/env python3
from meta_analysis import flip_strand, is_symmetric, format_num
import datetime
import argparse
import json
import gzip
from collections import namedtuple, defaultdict
import sys
import math
import scipy.stats
import numpy
from typing import Dict, Tuple, List
import subprocess
from collections import deque
import re

flip = {"A":"T","C":"G","T":"A","G":"C"}


class Variant():

    def __init__(self, chr, pos, ref, alt, af, filt, an):
        self.chr = int(chr)
        self.pos = int(float(pos))
        self.ref = ref.strip().upper()
        self.alt = alt.strip().upper()
        self.af = float(af) if af != "NA" else None
        self.filt = filt
        self.an = int(an)

    def __eq__(self, other):

        return self.chr == other.chr and self.pos == other.pos and self.ref == other.ref and self.alt == other.alt

    def __lt__(self, other):

        return (  (self.chr==other.chr and self.pos<other.pos)
                  or (self.chr < other.chr)
               )

    def is_equal(self, other:'Variant') -> bool:
        """
            Checks if this Variant is the same variant (possibly different strand or ordering of alleles)
            returns: true if the same false if not
        """

        if (self.chr == other.chr and self.pos == other.pos):
            flip_ref =  flip_strand(other.ref)
            flip_alt =  flip_strand(other.alt)

            if self.ref== other.ref and self.alt == other.alt :
                return True

            if is_symmetric( other.ref, other.alt ):
                ## never strandflip symmetrics. Assumed to be aligned.
                if self.ref == other.ref and self.alt == other.alt:
                    return True
                elif self.ref == other.alt and self.alt == other.ref:
                    return True

            elif (self.ref == other.alt and self.alt == other.ref) :
                return True
            elif (self.ref == flip_ref and self.alt==flip_alt):
                return True
            elif (self.ref == flip_alt and self.alt==flip_ref):
                return True

        return False
    

class VariantData(Variant):

    def __init__(self, chr, pos, ref, alt, af, beta, extra_cols=[], gnomad_af=None, af_fc=None, gnomad_filt="NA"):
        self.chr = int(chr)
        self.pos = int(float(pos))
        self.ref = ref.strip().upper()
        self.alt = alt.strip().upper()
        self.af = float(af) if af != 'NA' else None
        self.beta = float(beta) if beta != 'NA' else None
        self.extra_cols = extra_cols
        self.gnomad_af = gnomad_af
        self.af_fc = af_fc
        self.gnomad_filt = gnomad_filt

    def equalize_to(self, other:'Variant') -> bool:
        """
            Checks if this VariantData is the same variant as given other variant (possibly different strand or ordering of alleles)
            If it is, changes this variant's alleles, beta and af accordingly
            returns: true if the same (flips effect direction, ref/alt alleles and af if necessary) or false if not the same variant
        """

        if (self.chr == other.chr and self.pos == other.pos):
            flip_ref =  flip_strand(other.ref)
            flip_alt =  flip_strand(other.alt)

            if self.ref== other.ref and self.alt == other.alt :
                return True

            if is_symmetric( other.ref, other.alt ):
                ## never strandflip symmetrics. Assumed to be aligned.
                if self.ref == other.ref and self.alt == other.alt:
                    return True
                elif self.ref == other.alt and self.alt == other.ref:
                    self.beta = -1 * self.beta if self.beta is not None else None
                    self.af = 1 - self.af if self.af is not None else None
                    t = self.alt
                    self.alt = self.ref
                    self.ref = t
                    return True

            elif (self.ref == other.alt and self.alt == other.ref) :
                self.beta = -1 * self.beta if self.beta is not None else None
                self.af = 1 - self.af if self.af is not None else None
                t = self.alt
                self.alt = self.ref
                self.ref = t
                return True
            elif (self.ref == flip_ref and self.alt==flip_alt):
                self.ref = flip_strand(self.ref)
                self.alt = flip_strand(self.alt)
                return True
            elif (self.ref == flip_alt and self.alt==flip_ref):
                self.beta = -1 * self.beta if self.beta is not None else None
                self.af = 1 - self.af if self.af is not None else None
                self.ref = flip_strand(self.alt)
                self.alt = flip_strand(self.ref)
                return True

        return False

    def __str__(self):
        cols = [str(self.chr), str(self.pos), self.ref, self.alt, str(self.af), str(self.beta)]
        cols.extend(self.extra_cols)
        cols.extend([format_num(self.gnomad_af, 3), format_num(self.af_fc, 3), self.gnomad_filt])
        return '\t'.join(cols)

def harmonize(file_in, file_ref, chr_col, pos_col, ref_col, alt_col, af_col, beta_col, require_gnomad, passing_only, gnomad_min_an, gnomad_max_abs_diff):

    required_cols = [chr_col, pos_col, ref_col, alt_col, af_col, beta_col]
    
    fp_ref = gzip.open(file_ref, 'rt')
    ref_has_lines = True
    ref_chr = 1
    ref_pos = 0
    ref_h_idx = {h:i for i,h in enumerate(fp_ref.readline().strip().split('\t'))}

    with gzip.open(file_in, 'rt') as f:
        header = f.readline().strip().split('\t')
        h_idx = {h:i for i,h in enumerate(header)}
        extra_cols = [h for h in header if not h in required_cols]
        print('\t'.join(required_cols + extra_cols + ['af_gnomad','af_fc','filt_gnomad']))
        for line in f:
            s = line.strip().split('\t')
            var = VariantData(chr = s[h_idx[chr_col]].replace('chr', '').replace('X', '23').replace('Y', '24'),
                              pos = s[h_idx[pos_col]],
                              ref = s[h_idx[ref_col]],
                              alt = s[h_idx[alt_col]],
                              af = s[h_idx[af_col]],
                              beta = s[h_idx[beta_col]],
                              extra_cols = [s[h_idx[extra_col]] for extra_col in extra_cols])
            ref_vars = []
            while ref_has_lines and ref_chr < var.chr or (ref_chr == var.chr and ref_pos < var.pos):
                ref_line = fp_ref.readline()
                if ref_line != '':
                    r = ref_line.strip().split('\t')
                    ref_chr = int(ref_line[ref_h_idx['#chr']])
                    ref_pos = int(ref_line[ref_h_idx['pos']])
                else:
                    ref_has_lines = False

            while ref_has_lines and ref_chr == var.chr and ref_pos == var.pos:
                ref_vars.append(Variant(chr = ref_chr,
                                        pos = ref_pos,
                                        ref = ref_line[ref_h_idx['ref']],
                                        alt = ref_line[ref_h_idx['alt']],
                                        af = ref_line[ref_h_idx['af_alt']],
                                        filt = ref_line[ref_h_idx['filter']],
                                        an = ref_line[ref_h_idx['an']]))
                ref_line = fp_ref.readline()
                if ref_line != '':
                    r = ref_line.strip().split('\t')
                    ref_chr = int(ref_line[ref_h_idx['#chr']])
                    ref_pos = int(ref_line[ref_h_idx['pos']])
                else:
                    ref_has_lines = False

            equal = []
            diffs = []
            fcs = []
            for r in ref_vars:
                if var == r and (not passing_only or r.filt == 'PASS') and r.an >= gnomad_min_an:
                    diff = 1
                    fc = 1e9
                    if r.af is not None and var.af is not None:
                        diff = abs(var.af - float(r.af))
                        fc = var.af/float(r.af) if float(r.af) != 0 else 1e9
                    equal.append(r)
                    diffs.append(diff)
                    fcs.append(fc)

            best_diff = -1
            if len(equal) > 0:
                best_diff = 2
                for i,diff in enumerate(diffs):
                    if diff < best_diff or (diff == best_diff and equal[i].ref == var.ref and equal[i].alt == var.alt):
                        best_diff = diff
                        best_diff_idx = i
                var.equalize_to(equal[best_diff_idx])
                var.gnomad_af = equal[best_diff_idx].af
                var.af_fc = fcs[best_diff_idx] if fcs[best_diff_idx] != 1e9 else None
                var.gnomad_filt = equal[best_diff_idx].filt

            if (not require_gnomad or len(equal) > 0) and best_diff <= gnomad_max_abs_diff:
                print(var)
            
def run():
    parser = argparse.ArgumentParser(description="Harmonize GWAS summary stats to reference")
    parser.add_argument('file_in', action='store', type=str, help='GWAS summary stats')
    parser.add_argument('file_ref', action='store', type=str, help='gnomAD reference file')
    parser.add_argument('--chr_col', action='store', type=str, default='#CHR', help='Chromosome column')
    parser.add_argument('--pos_col', action='store', type=str, default='POS', help='Position column')
    parser.add_argument('--ref_col', action='store', type=str, default='REF', help='Reference allele column')
    parser.add_argument('--alt_col', action='store', type=str, default='ALT', help='Alternative allele column')
    parser.add_argument('--af_col', action='store', type=str, default='af_alt', help='Allele frequency allele column')
    parser.add_argument('--beta_col', action='store', type=str, default='beta', help='Beta column')
    parser.add_argument('--require_gnomad', action='store_true', help='Filter out variants not in gnomAD')
    parser.add_argument('--passing_only', action='store_true', help='Filter out non-passing variants in gnomAD')
    parser.add_argument('--gnomad_min_an', action='store', type=int, default=0, help='Minimum AN in gnomAD')
    parser.add_argument('--gnomad_max_abs_diff', action='store', type=float, default=1.0, help='Maximum absolute difference between variant and gnomAD AF')
    args = parser.parse_args()
    harmonize(file_in = args.file_in,
              file_ref = args.file_ref,
              chr_col = args.chr_col,
              pos_col = args.pos_col,
              ref_col = args.ref_col,
              alt_col = args.alt_col,
              af_col = args.af_col,
              beta_col = args.beta_col,
              require_gnomad = args.require_gnomad,
              passing_only = args.passing_only,
              gnomad_min_an = args.gnomad_min_an,
              gnomad_max_abs_diff = args.gnomad_max_abs_diff)
    
if __name__ == '__main__':
    run()