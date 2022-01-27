#!/usr/bin/env python3
from meta_analysis import flip_strand, is_symmetric, format_num
import argparse
import gzip
import numpy
from typing import Dict, Tuple, List

flip = {'A':'T','C':'G','T':'A','G':'C'}


class Variant():

    def __init__(self, chr, pos, ref, alt, af):
        self.chr = int(str(chr).replace('chr', '').replace('X', '23').replace('Y', '24'))
        self.pos = int(float(pos))
        self.ref = ref.strip().upper()
        self.alt = alt.strip().upper()
        self.af = float(af) if af != "NA" else None

    def __eq__(self, other):
        return self.chr == other.chr and self.pos == other.pos and self.ref == other.ref and self.alt == other.alt

    def __lt__(self, other):
        return self.chr < other.chr or (self.chr == other.chr and self.pos < other.pos)

    def is_equal(self, other:'Variant') -> bool:
        """Checks if two variants are equal.

        Checks if this Variant is the same variant as given other variant (possibly different strand or ordering of alleles).
        
        Args:
            other:
                Variant as Variant object to compare this variant to.

        Returns:
            True if the same or False if not the same variant.
        """

        if self.chr == other.chr and self.pos == other.pos:

            if self.ref == other.ref and self.alt == other.alt :
                return True

            flip_ref = flip_strand(other.ref)
            flip_alt = flip_strand(other.alt)

            if is_symmetric(other.ref, other.alt):
                ## never strandflip symmetrics. Assumed to be aligned.
                return False
            elif self.ref == other.alt and self.alt == other.ref:
                return True
            elif self.ref == flip_ref and self.alt == flip_alt:
                return True
            elif self.ref == flip_alt and self.alt == flip_ref:
                return True

        return False
    

class VariantData(Variant):

    def __init__(self, chr, pos, ref, alt, af, beta, extra_cols=[], gnomad_af=None, gnomad_filt='NA'):
        super().__init__(chr, pos, ref, alt, af)
        self.beta = float(beta) if beta != 'NA' else None
        self.extra_cols = extra_cols
        self.gnomad_af = gnomad_af
        self.fc = None
        self.gnomad_filt = gnomad_filt
        self.af_diff = None

    def equalize_to(self, other:'Variant') -> bool:
        """Checks if two variants are equal and changes this variant's alleles and beta accordingly.

        Checks if this Variant is the same variant as given other variant (possibly different strand or ordering of alleles) and changes this variant's alleles and beta accordingly.
        
        Args:
            other:
                Variant as Variant object to compare this variant to.

        Returns:
            True if the same or False if not the same variant.
        """

        if self.chr == other.chr and self.pos == other.pos:

            if self.ref== other.ref and self.alt == other.alt :
                return True
            
            flip_ref = flip_strand(other.ref)
            flip_alt = flip_strand(other.alt)

            if is_symmetric( other.ref, other.alt ):
                ## never strandflip symmetrics. Assumed to be aligned.
                return False

            elif self.ref == other.alt and self.alt == other.ref:
                self.beta = -1 * self.beta if self.beta is not None else None
                self.af = 1 - self.af if self.af is not None else None
                t = self.alt
                self.alt = self.ref
                self.ref = t
                return True
            elif self.ref == flip_ref and self.alt==flip_alt:
                self.ref = flip_strand(self.ref)
                self.alt = flip_strand(self.alt)
                return True
            elif self.ref == flip_alt and self.alt==flip_ref:
                self.beta = -1 * self.beta if self.beta is not None else None
                self.af = 1 - self.af if self.af is not None else None
                self.ref = flip_strand(self.alt)
                self.alt = flip_strand(self.ref)
                return True

        return False

    @property
    def af_fc(self):
        if self.fc is None and self.gnomad_af is not None and self.af is not None:
            self.fc = self.af/self.gnomad_af if self.gnomad_af != 0 else 1e9
        return self.fc

    @property
    def af_abs_diff(self):
        if self.af_diff is None and self.gnomad_af is not None and self.af is not None:
            self.af_diff = abs(self.af-self.gnomad_af)
        elif self.af_diff is None and (self.gnomad_af is None or self.af is None):
            self.af_diff = -1
        return self.af_diff

    def __str__(self):
        cols = [str(self.chr), str(self.pos), self.ref, self.alt, str(self.af), str(self.beta)]
        cols.extend(self.extra_cols)
        cols.extend([format_num(self.gnomad_af, 3), format_num(self.af_fc, 3), self.gnomad_filt])
        return '\t'.join(cols)


class VariantGnomad(Variant):

    def __init__(self, chr, pos, ref, alt, af, filt, an):
        super().__init__(chr, pos, ref, alt, af)
        self.filt = filt
        self.an = int(an)


def choose_best_variant(variant_list: List[VariantData], require_gnomad: bool, gnomad_max_abs_diff: float) -> VariantData:
    """Finds best match according to AF difference from list of variants

    Args:
        variant_list:
            List of VariantData objects
        require_gnomad:
            If variant not found in gnomAD, don't print variant
        gnomad_max_abs_diff:
            Maximum absolute difference between variant and gnomAD AF

    Returns:
        VarianData object containing the best match or None if variant does not pass filters
    """

    if len(variant_list) > 1:
        diffs = [var.af_abs_diff for var in variant_list]
        best_idx = diffs.index(min(diffs))
        if (not require_gnomad or variant_list[best_idx].gnomad_af is not None) and variant_list[best_idx].af_abs_diff <= gnomad_max_abs_diff:
            return variant_list[best_idx]
    elif (not require_gnomad or variant_list[0].gnomad_af is not None) and variant_list[0].af_abs_diff <= gnomad_max_abs_diff:
        return variant_list[0]
    else:
        return None


def harmonize(file_in, file_ref, chr_col, pos_col, ref_col, alt_col, af_col, beta_col, require_gnomad, passing_only, gnomad_min_an, gnomad_max_abs_diff, pre_aligned, keep_best_duplicate):

    required_cols = [chr_col, pos_col, ref_col, alt_col, af_col, beta_col]
    
    fp_ref = gzip.open(file_ref, 'rt')
    ref_has_lines = True
    ref_chr = 1
    ref_pos = 0
    ref_h_idx = {h:i for i,h in enumerate(fp_ref.readline().strip().split('\t'))}

    previous_vars = []
    ref_vars = []

    with gzip.open(file_in, 'rt') as f:
        header = f.readline().strip().split('\t')
        h_idx = {h:i for i,h in enumerate(header)}
        extra_cols = [h for h in header if not h in required_cols]
        print('\t'.join(required_cols + extra_cols + ['af_gnomad','af_fc','filt_gnomad']))
        for line in f:
            s = line.strip().split('\t')
            var = VariantData(chr = s[h_idx[chr_col]],
                              pos = s[h_idx[pos_col]],
                              ref = s[h_idx[ref_col]],
                              alt = s[h_idx[alt_col]],
                              af = s[h_idx[af_col]],
                              beta = s[h_idx[beta_col]],
                              extra_cols = [s[h_idx[extra_col]] for extra_col in extra_cols])
            
            if not ref_vars or ref_vars[0] < var:
                ref_vars = []
                while ref_has_lines and (ref_chr < var.chr or (ref_chr == var.chr and ref_pos < var.pos)):
                    ref_line = fp_ref.readline()
                    if ref_line != '':
                        r = ref_line.strip().split('\t')
                        ref_chr = int(r[ref_h_idx['#chr']])
                        ref_pos = int(r[ref_h_idx['pos']])
                    else:
                        ref_has_lines = False

                while ref_has_lines and ref_chr == var.chr and ref_pos == var.pos:
                    ref_vars.append(VariantGnomad(chr = ref_chr,
                                                pos = ref_pos,
                                                ref = r[ref_h_idx['ref']],
                                                alt = r[ref_h_idx['alt']],
                                                af = r[ref_h_idx['af_alt']],
                                                filt = r[ref_h_idx['filter']],
                                                an = r[ref_h_idx['an']]))
                    ref_line = fp_ref.readline()
                    if ref_line != '':
                        r = ref_line.strip().split('\t')
                        ref_chr = int(r[ref_h_idx['#chr']])
                        ref_pos = int(r[ref_h_idx['pos']])
                    else:
                        ref_has_lines = False

            equal = []
            diffs = []
            for ref_var in ref_vars:
                if (var == ref_var or (not pre_aligned and var.is_equal(ref_var))) and (not passing_only or ref_var.filt == 'PASS') and ref_var.an >= gnomad_min_an:
                    diff = 1
                    if ref_var.af is not None and var.af is not None:
                        if not pre_aligned:
                            var.equalize_to(ref_var)
                        diff = abs(var.af - ref_var.af)
                    equal.append(ref_var)
                    diffs.append(diff)

            best_diff = -1
            if len(equal) > 0:
                best_diff = 2
                for i,diff in enumerate(diffs):
                    if diff < best_diff or (diff == best_diff and equal[i].ref == var.ref and equal[i].alt == var.alt):
                        best_diff = diff
                        best_diff_idx = i
                if not pre_aligned:
                    var.equalize_to(equal[best_diff_idx])
                var.gnomad_af = equal[best_diff_idx].af
                var.gnomad_filt = equal[best_diff_idx].filt

            if previous_vars and previous_vars[0] != var:
                if len(previous_vars) == 1 or keep_best_duplicate:
                    best_var = choose_best_variant(variant_list=previous_vars, require_gnomad=require_gnomad, gnomad_max_abs_diff=gnomad_max_abs_diff)
                    if best_var:
                        print(best_var)
                        best_var = None
                previous_vars = []
            previous_vars.append(var)

        if previous_vars:
            if len(previous_vars) == 1 or keep_best_duplicate:
                best_var = choose_best_variant(variant_list=previous_vars, require_gnomad=require_gnomad, gnomad_max_abs_diff=gnomad_max_abs_diff)
                if best_var:
                    print(best_var)
                    best_var = None


def run():
    parser = argparse.ArgumentParser(description="Harmonize GWAS summary stats to reference")
    parser.add_argument('file_in', action='store', type=str, help='GWAS summary stats')
    parser.add_argument('file_ref', action='store', type=str, help='GnomAD reference file')
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
    parser.add_argument('--pre_aligned', action='store_true', help='Input summary stats are already aligned to reference (disables flipping of alleles to try find best match)')
    parser.add_argument('--keep_best_duplicate', action='store_true', help='If duplicate variants (by chr:pos:ref:alt) keep variant with smallest AF difference to reference. Otherwise discard both')
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
              gnomad_max_abs_diff = args.gnomad_max_abs_diff,
              pre_aligned = args.pre_aligned,
              keep_best_duplicate = args.keep_best_duplicate)
    
if __name__ == '__main__':
    run()