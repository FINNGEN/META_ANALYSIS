#!/usr/bin/env python3
import argparse
import json
import gzip
import sys
import math
from scipy.stats import chi2
import scipy.stats
import numpy
from typing import Dict, Tuple, List
import subprocess
from collections import deque
import re

chrord = {"chr"+str(chr): chr for chr in range(1, 23)}
chrord.update({"X": 23, "Y": 24, "MT": 25, "chrX": 23, "chrY": 24, "chrMT": 25})
chrord.update({str(chr): chr for chr in range(1, 25)})
chrord.update({chr: chr for chr in range(1, 25)})

re_allele = re.compile('^[ATCG]+$', re.IGNORECASE)

def het_test( effs_sizes: List[float], weights: List[float], effs_size_meta: float) -> float:
    '''
        Computes Cochran's Q test for heterogeneity
        input:
            effs_sizes: original effect sizes
            weights: weights
            effs_size_meta: effect size from meta-analysis
        output:
            p-value
    '''
    k=len(effs_sizes)

    effs_sizes_array=numpy.array(effs_sizes)
    weights_array=numpy.array(weights)    
    eff_dev=weights_array*((effs_sizes_array-effs_size_meta)**2)
    sum_eff_dev=numpy.sum(eff_dev)

    return scipy.stats.distributions.chi2.sf(sum_eff_dev, k-1)

def n_meta( studies : List[Tuple['Study','VariantData']] ) -> Tuple:
    '''
        Computes sample size weighted meta-analysis for variants in studies
        input:
            studies: studies and data in tuples
        output:
            tuple with results from meta-analysis or None
    '''
    weights = []
    effs_size_org = [] 

    effs_size = []
    tot_size =0
    sum_betas=0
    sum_weights=0
    for s in studies:
        study = s[0]
        dat = s[1]
        effs_size.append( math.sqrt(study.effective_size) * numpy.sign(dat.beta) * dat.z_score)
        sum_weights+= math.sqrt(study.effective_size)
        sum_betas+=math.sqrt(study.effective_size) * dat.beta
        tot_size+=study.effective_size
        weights.append(math.sqrt(study.effective_size))
        effs_size_org.append(dat.beta)

    beta_meta=sum_betas/sum_weights

    z_val_abs = abs(sum(effs_size)) / math.sqrt(tot_size)
    
    #TODO se
    return (
        beta_meta,
        None,
        max(sys.float_info.min * sys.float_info.epsilon, 2 * scipy.stats.norm.sf(z_val_abs)),
        -(scipy.stats.norm.logsf(z_val_abs) + math.log(2)) / math.log(10),
        effs_size_org,
        weights
    ) if len(effs_size)==len(studies) else None


def inv_var_meta( studies : List[Tuple['Study','VariantData']] ) -> Tuple:
    '''
        Computes inverse-variance weighted meta-analysis for variants in studies
        input:
            studies: studies and data in tuples
        output:
            tuple with results from meta-analysis or None
    '''
    weights = []
    effs_size_org = []

    effs_inv_var = []
    sum_inv_var=0
    for s in studies:
        study = s[0]
        dat = s[1]
        if dat.se is None or dat.se==0:
            print("Standard error was none/zero for variant " + str(dat) + " in study " + study.name, file=sys.stderr)
            return None
        var = (dat.se * dat.se)

        inv_var =  (1/var)
        sum_inv_var+=inv_var
        effs_inv_var.append( inv_var *  dat.beta )

        weights.append(inv_var)
        effs_size_org.append(dat.beta)

    beta_meta=sum(effs_inv_var)/ sum_inv_var

    z_val_abs = abs(sum(effs_inv_var) / math.sqrt(sum_inv_var))
    
    return (
        beta_meta,
        math.sqrt(1/sum_inv_var),
        max(sys.float_info.min * sys.float_info.epsilon, 2 * scipy.stats.norm.sf(z_val_abs)),
        -(scipy.stats.norm.logsf(z_val_abs) + math.log(2)) / math.log(10),
        effs_size_org,
        weights
    )


def variance_weight_meta( studies : List[Tuple['Study','VariantData']] ) -> Tuple:
    '''
        Computes variance weighted meta-analysis for variants in studies
        input:
            studies: studies and data in tuples
        output:
            tuple with results from meta-analysis or None
    '''
    weights = []
    effs_size_org = []

    effs_se = []
    tot_se = 0
    sum_weights=0
    sum_betas=0
    for s in studies:
        study = s[0]
        dat = s[1]

        if dat.se is None or dat.se==0:
            print("Standard error was none/zero for variant " + str(dat) + " in study " + study.name, file=sys.stderr)
            return None
        weight =  (1/dat.se) * dat.z_score
        sum_weights+=weight
        sum_betas+= weight * dat.beta
        effs_se.append( weight * numpy.sign(dat.beta)  )
        tot_se+=1/ (dat.se * dat.se)
        
        weights.append(weight)
        effs_size_org.append(dat.beta)

    beta_meta=sum_betas / sum_weights

    z_val_abs = abs(sum(effs_se)) / math.sqrt(tot_se)
    
    #TODO SE
    return (
        beta_meta,
        None,
        max(sys.float_info.min * sys.float_info.epsilon, 2 * scipy.stats.norm.sf(z_val_abs)),
        -(scipy.stats.norm.logsf(z_val_abs) + math.log(2)) / math.log(10),
        effs_size_org,
        weights
    )

SUPPORTED_METHODS = {"n":n_meta,"inv_var":inv_var_meta,"variance":variance_weight_meta}

def check_eff_field(field):
    if field.lower() in ["beta","or"]:
        return field.lower()
    else:
        raise Exception("effect_type must be beta or OR")

flip = {"A":"T","C":"G","T":"A","G":"C"}

def flip_strand( allele):
    return "".join([ flip[a] for a in allele])

def is_symmetric(a1, a2):
    return (a1=="A" and a2=="T") or (a1=="T" and a2=="A") or (a1=="C" and a2=="G") or (a1=="G" and a2=="C")


class VariantData:

    def __init__(self, chr, pos, ref, alt, beta, pval, se=None, extra_cols=[]):
        self.chr = chr
        self.pos = int(float(pos))
        self.ref = ref.strip().upper()
        self.alt = alt.strip().upper()
        self.beta = beta
        self.pval = pval
        self.z_scr = None
        self.indel = None
        try:
            self.se = float(se) if se is not None  else None
        except ValueError:
            self.se = None

        self.extra_cols = extra_cols

    def __eq__(self, other):

        return self.chr == other.chr and self.pos == other.pos and self.ref == other.ref and self.alt == other.alt

    def __lt__(self, other):

        return (  (self.chr==other.chr and self.pos<other.pos)
                  or (self.chr < other.chr)
               )

    def _equalizer(self, other: 'VariantData', equalize: bool = False, flip_indels: bool = False) -> bool:
        """Checks if two variants are equal

        Checks if this VariantData is the same variant as given other variant (possibly different strand or ordering of alleles).
        
        Args:
            other:
                Variant as VariantData object to compare this variant to.
            equalize:
                If True changes this variant's alleles and beta accordingly.
            flip_indels:
                If True will try matching indels with flipping.

        Returns:
            True if the same or False if not the same variant.
        """

        if (self.chr == other.chr and self.pos == other.pos):

            if self.ref == other.ref and self.alt == other.alt :
                return True
            
            if not flip_indels and self.is_indel:
                return False
            
            flip_ref = flip_strand(other.ref)
            flip_alt = flip_strand(other.alt)

            if is_symmetric(other.ref, other.alt):
                ## never strandflip symmetrics. Assumed to be aligned.
                return False

            elif (self.ref == other.alt and self.alt == other.ref) :
                if equalize:
                    self.beta = -1 * self.beta if self.beta is not None else None
                    t = self.alt
                    self.alt = self.ref
                    self.ref = t
                return True
            elif (self.ref == flip_ref and self.alt == flip_alt):
                if equalize:
                    self.ref = flip_strand(self.ref)
                    self.alt = flip_strand(self.alt)
                return True
            elif (self.ref == flip_alt and self.alt == flip_ref):
                if equalize:
                    self.beta = -1 * self.beta if self.beta is not None else None
                    self.ref = flip_strand(self.alt)
                    self.alt = flip_strand(self.ref)
                return True

        return False

    def is_equal(self, other: 'VariantData', flip_indels: bool = False) -> bool:
        """Checks if two variants are equal

        Checks if this VariantData is the same variant as given other variant (possibly different strand or ordering of alleles).
        
        Args:
            other:
                Variant as VariantData object to compare this variant to.
            flip_indels:
                If True will try matching indels with flipping.

        Returns:
            True if the same or False if not the same variant.
        """
        return self._equalizer(other=other, equalize=False, flip_indels=flip_indels)

    def equalize_to(self, other: 'VariantData', flip_indels: bool = False) -> bool:
        """Checks if two variants are equal and changes this variant's alleles and beta accordingly

        Checks if this VariantData is the same variant as given other variant (possibly different strand or ordering of alleles) and changes this variant's alleles and beta accordingly.
        
        Args:
            other:
                Variant as VariantData object to compare this variant to.
            flip_indels:
                If True will try matching indels with flipping.

        Returns:
            True if the same or False if not the same variant.
        """
        return self._equalizer(other=other, equalize=True, flip_indels=flip_indels)

    @property
    def z_score(self):
        '''
            Lazy compute unsigned z-score
        '''
        if self.z_scr is None:
            self.z_scr = math.sqrt(chi2.isf(self.pval, df=1))
        return self.z_scr

    @property
    def is_indel(self):
        if self.indel is None:
            self.indel = len(self.ref)>1 or len(self.alt)>1
        return self.indel

    def __str__(self):
        return "chr:{} pos:{} ref:{} alt:{} beta:{} pval:{} se:{} ".format(self.chr, self.pos, self.ref, self.alt, self.beta, self.pval, self.se)


class Study:
    REQUIRED_DATA_FIELDS = {"chr":str,"pos":str,"ref":str,"alt":str, "effect":str,
    "pval":str}

    REQUIRED_CONF = {"name":str,"file":str, "n_cases": int, "n_controls":int,
    "chr":str,"pos":str,"ref":str,"alt":str, "effect":str,
    "effect_type":check_eff_field,
    "pval":str}

    OPTIONAL_FIELDS = {"se":str}

    def __init__(self, conf, chrom=None, sep='\t', flip_indels=False):
        '''
        chrom: a chromosome to limit to or None if all chromosomes
        sep: field separator (default: "\t")
        '''
        self.conf = conf
        self.chrom = chrord[chrom] if chrom is not None else None
        self.sep = sep
        self.flip_indels = flip_indels
        self.future = deque()
        self.eff_size = None
        self.z_scr = None
        self.prev_var = None
        self.header = None
        for v in Study.REQUIRED_CONF:
            if v not in self.conf:
                raise Exception("Meta configuration for study must contain required elements: "
                    + ",".join(Study.REQUIRED_CONF.keys() ) + ". Offending configuration: " + str(self.conf))

            try:
                self.conf[v] = Study.REQUIRED_CONF[v](self.conf[v])
            except Exception as e:
                raise Exception("Illegal data type in configuration for field " + str(v) +
                    " in configuration: " + str(self.conf) + ". ERR:" + str(e))

        for v in Study.OPTIONAL_FIELDS:
            if v not in self.conf:
                continue
            try:
                self.conf[v] = Study.OPTIONAL_FIELDS[v](self.conf[v])
            except Exception as e:
                raise Exception("Illegal data type in configuration for field " + v +
                    " in configuration: " + str(self.conf) + ". ERR:" + str(e))

        self.conf["fpoint"] = gzip.open(conf["file"],'rt')
        self.header = conf["fpoint"].readline().rstrip().split(self.sep)

        for k in Study.REQUIRED_DATA_FIELDS.keys():
            if self.conf[k] not in self.header:
                raise Exception("Required headers not in data in study " + self.conf["name"] + ". Missing:" + ",".join([ self.conf[k] for k in Study.REQUIRED_DATA_FIELDS.keys() if self.conf[k] not in self.header])  )
        self.conf["h_idx"] = { k:self.header.index( self.conf[k] ) for k in Study.REQUIRED_DATA_FIELDS.keys() }

        for f in Study.OPTIONAL_FIELDS.keys():
            if f in self.conf:
                 if self.conf[f] not in self.header:
                     raise Exception("Configured column " + self.conf[f] + " not found in the study results " + self.conf["name"])
                 self.conf["h_idx"][f] = self.header.index(self.conf[f])

        if "extra_cols" in self.conf:
            for c in self.conf["extra_cols"]:
                if c not in self.header:
                    raise Exception("Configured column " + self.conf[c] + " not found in the study results " + self.conf["name"])
                self.conf["h_idx"][c] = self.header.index(c)
        else:
             self.conf["extra_cols"] = []



    @property
    def n_cases(self):
        return self.conf["n_cases"]

    @property
    def n_controls(self):
        return self.conf["n_controls"]

    @property
    def effective_size(self):
        if self.eff_size is None:
            self.eff_size = ( (4 * self.n_cases *  self.n_controls  ) / ( self.n_cases+  self.n_controls ))
        return self.eff_size

    @property
    def name(self):
        return self.conf["name"]

    def has_std_err(self):
        return "se" in self.conf

    def get_next_data(self, just_one: bool = False) -> List[VariantData]:
        """
            Returns a list of variants. List containts >1 elements if they are on the same position and just_one == False.

            input:
                just_one: always returns only the next variant in order and not all next with the same position
            output:
                list of next variants
        """

        vars = list()

        if len(self.future)>0:
            ## only return variants with same position so that possible next variant position stored stays
            f = [ (i,v) for i,v in enumerate(self.future) if i==0 or (v.chr==self.future[0].chr and v.pos==self.future[0].pos) ]
            for i,v in reversed(f):
                del self.future[i]
            vars.extend([ v for i,v in f ])
            if len(self.future)>0:
                return vars

        while True:
            chr = None
            ref = None
            alt = None
            l = None
            ## loop ignoring  alternate contigs and non-ATCG alleles for now.
            while chr is None or chr not in chrord or (self.chrom is not None and chr != self.chrom) or re_allele.match(ref) is None or re_allele.match(alt) is None:
                l = self.conf["fpoint"].readline()
                if l=="":
                    return None if len(vars) == 0 else vars

                l = l.rstrip().split(self.sep)

                if len(l) != len(self.header):
                    raise Exception("Number of fields in header and line do not match for study " + self.name + " in file " + self.conf['file'] + ".\nOffending line: " + "\t".join(l))

                chr = l[self.conf["h_idx"]["chr"]]
                chr = chrord[chr] if chr in chrord else None
                ref = l[self.conf["h_idx"]["ref"]]
                alt = l[self.conf["h_idx"]["alt"]]

            pos = l[self.conf["h_idx"]["pos"]]
            eff = l[self.conf["h_idx"]["effect"]]
            pval = l[self.conf["h_idx"]["pval"]]

            se = l[self.conf["h_idx"]["se"]] if "se" in self.conf["h_idx"] else None

            effect_type = self.conf["effect_type"]
            try:
                pval = float(pval)
                eff = float(eff)
            except:
                pval = None
                eff = None

            if( effect_type=="or" and eff):
                eff = math.log(eff)

            extracols = [ l[self.conf["h_idx"][c]] for c in self.conf["extra_cols"] ]

            v = VariantData(chr,pos,ref,alt, eff, pval, se, extracols)

            if self.prev_var is not None and v < self.prev_var:
                raise Exception("Disorder in study " + self.name + " in file " + self.conf['file'] + ". Sort all summary statistic files by chromosome and then position and rerun.\nOffending line: " + "\t".join(l))
            self.prev_var = v
            if len(vars)==0 or ( vars[0].chr == v.chr and vars[0].pos == v.pos  ):
                added=False
                for v_ in vars:
                    if v == v_:
                        print('ALREADY ADDED FOR STUDY ' + self.name + ': ' + str(v), file=sys.stderr)
                        added=True
                if not added:
                    vars.append(v )
                if just_one:
                    break
            else:
                self.future.append(v )
                break

        return vars

    @property
    def extra_cols(self):
        return self.conf["extra_cols"]

    def get_match(self, dat: VariantData) -> VariantData:
        """
            Reads current study until variant in 'dat' is reached or overtaken in chr pos orded.
            IF matching variant found (can flip alleles) the matching VariantData(effect flipped if alleles flipped) is returned.
            input:
                dat: the variant to look for
            output: matching VariantData in this study or None if no match.
        """

        otherdats = self.get_next_data( )

        if otherdats is None or len(otherdats)==0:
            return None

        while otherdats is not None and (otherdats[0].chr<dat.chr or (otherdats[0].chr==dat.chr and otherdats[0].pos<dat.pos)):
            otherdats = self.get_next_data()

        if otherdats is None:
            return None

        if otherdats[0].chr > dat.chr or otherdats[0].pos> dat.pos:
            self.put_back(otherdats)
            return None

        for i,v in enumerate(otherdats):
            if v.equalize_to(dat, flip_indels=self.flip_indels):
                del otherdats[i]
                self.put_back(otherdats)
                return v

        ## no match but stayed in the same pos. add variants back to future queue
        self.put_back(otherdats)
        return None


    def put_back(self, variantlist: List[VariantData]):
        '''
        Put list of variants back to wait for matching
        
        input:
            variantlist: list of VariantData objects
        output:
            p-value
        '''
        
        self.future.extendleft(variantlist)


def get_studies(conf:str, chrom, sep, flip_indels) -> List[Study]:
    """
        Reads json configuration and returns studies in the meta
    """

    studies_conf = json.load(open(conf,'r'))

    return [ Study(s, chrom, sep, flip_indels) for s in studies_conf["meta"]]

def format_num(num, precision=2):
    return "NA" if num is None or numpy.isnan(num) else numpy.format_float_scientific(num, precision=precision)

def do_meta(study_list: List[ Tuple[Study, VariantData]], methods: List[str], is_het_test) -> List[Tuple] :
    '''
        Computes meta-analysis between all studies and data given in the std_list
        input:
            study_list: studies and data in tuples
            methods: list of methods to calculate
            is_het_test: boolean, do heterogeneity test
        output:
            list of tuples (effect_size, standard error, p-value (, het test p-value)) for each method in the same order as methods were given
    '''
    met = [ SUPPORTED_METHODS[m](study_list) for m in methods ]

    meta_res = []
    for m in met:
        if m is not None:
            if is_het_test:
                meta_res.append((format_num(m[0]), format_num(m[1]), format_num(m[2]), numpy.round(m[3], 2), format_num(het_test(m[4], m[5], m[0]))))
            else:
                meta_res.append((format_num(m[0]), format_num(m[1]), format_num(m[2]), numpy.round(m[3], 2)))
        else:
            meta_res.append(None)

    return meta_res

def get_next_variant( studies : List[Study]) -> List[VariantData]:
    '''
        get variant data for all studies
        The variant data is the first in chromosomal order across studies (ties broken by alphabetic order of ref)
        input:
            study_list: studies and data in tuples
        output:
            List of VariantData objects. The list is in the same order as the input studies. If smallest variant is
            not found in a study, that position in the list will be Null
    '''

    dats = []
    first = None
    for s in studies:
        d = s.get_next_data()
        dats.append(d)

        if d is not None:
            for v in d:
                if first is None or v < first:
                    first = v

    res = []
    for i,s in enumerate(studies):

        if dats[i] is None:
            res.append(None)
            continue

        # Flag tracks that only one variant (best match) is returned per study, if multiple variants equal to first
        added=False
        for j,v in reversed(list(enumerate(dats[i]))):
            if v == first:
                res.append(v)
                added=True
                del dats[i][j]
                s.put_back(dats[i])
                break
            if not v.is_equal(first, flip_indels=s.flip_indels):
                s.put_back([v])
                del dats[i][j]
        if not added:
            for j,v in reversed(list(enumerate(dats[i]))):
                if v.equalize_to(first, flip_indels=s.flip_indels):
                    res.append(v)
                    added=True
                    del dats[i][j]
                    break
            s.put_back(dats[i])
        if not added:
            res.append(None)

    return res

def validate_methods(methods, studies):
    M = []
    for m in methods:
        if m not in SUPPORTED_METHODS:
            raise Exception("Unsupported meta method" + m + " given. Supported values" + ",".join(SUPPORTED_METHODS))
        if m in ["inv_var", "variance"]:
            for s in studies:
                if not s.has_std_err():
                    raise Exception("Variance based method requested but not all studies have se column specified.")
        M.append(m)
    return M

def run():
    '''
        First parameter should be a path to a json configuration file with these elements:
            "name":"FINNGEN",
            "file":"/Users/mitja/projects/finngen/META_ANALYSIS/I9_AF.gz",
            "n_cases": 6570 ,
            "n_controls": 48378,
            "chr":"CHR",
            "pos":"POS",
            "ref":"Allele1",
            "alt":"Allele2",
            "effect":"BETA",
            "effect_type":"beta",
            "pval":"p.value"
            "se":"SE" <- this parameter is optional. If given for compared studies additional p-value will be added using this as a weight for z-score.
        Second parameter should be a path to (empty/not existing) directory where the data should be stored
    '''

    parser = argparse.ArgumentParser(description='Run x-way meta-analysis')
    parser.add_argument('config_file', action='store', type=str, help='Configuration file')
    parser.add_argument('path_to_res', action='store', type=str, help='Result file')
    parser.add_argument('methods', action='store', type=str, nargs="+", help='Methods to use in calculating meta-analysis statistics. Allowed values: n, inv_var, variance.', choices=["n", "inv_var", "variance"])

    parser.add_argument('--not_quiet', action='store_false', dest='quiet', help='Print matching variants to stdout')
    parser.add_argument('--leave_one_out', action='store_true', help='Do leave-one-out meta-analysis')
    parser.add_argument('--is_het_test', action='store_true', help='Do heterogeneity tests based on Cochrans Q and output het_p')
    parser.add_argument('--pairwise_with_first', action='store_true', help='Do pairwise meta-analysis with the first given study')
    parser.add_argument('--sep', default='\t', action='store', help='Input file field separator (Default: "\\t")')
    parser.add_argument('--chrom', action='store', type=str, help='Restrict to given chromosome')
    parser.add_argument('--flip_indels', action='store_true', help='Try variant aligning by flipping indels also. By default indels are not flipped')

    args = parser.parse_args()

    studs = get_studies(args.config_file, args.chrom, args.sep, args.flip_indels)
    methods = validate_methods(args.methods, studs)

    n_meta_cols = 5 if args.is_het_test else 4

    with open(args.path_to_res, 'w') as out:

        out.write("\t".join(["#CHR","POS","REF","ALT","SNP"]))

        ## align to leftmost STUDY
        for i,s in enumerate(studs):
            out.write( "\t" +  "\t".join( [ s.name + "_beta", s.name + "_sebeta", s.name + "_pval"] ))
            out.write( ("\t" if len(s.extra_cols) else "") + "\t".join( [s.name + "_" + c for c in s.extra_cols] ) )
            if args.pairwise_with_first and i>0:
                for m in methods:
                    out.write("\t" + studs[0].name + "_" + s.name + "_" +  m + "_meta_beta\t" + studs[0].name + "_" + s.name + "_" +  m + "_meta_sebeta\t" + studs[0].name + "_" + s.name + "_" +  m + "_meta_p\t" + studs[0].name + "_" + s.name + "_" +  m + "_meta_mlogp")

        out.write("\tall_meta_N")
        for m in methods:
            out.write("\tall_" + m + "_meta_beta\tall_" + m + "_meta_sebeta\tall_" + m + "_meta_p\tall_" + m + "_meta_mlogp")
            if args.is_het_test:
                out.write("\tall_" + m + "_het_p")

        if args.leave_one_out:
            for s in studs:
                out.write("\t" + "leave_" + s.name + "_N")
                for m in methods:
                    out.write( "\t" +  "\t".join( ["leave_" + s.name + "_" + m + "_meta_beta", "leave_" + s.name + "_" + m + "_meta_sebeta", "leave_" + s.name + "_" + m + "_meta_p", "leave_" + s.name + "_" + m + "_meta_mlogp"] ))
                    if args.is_het_test:
                        out.write("\tleave_" + s.name + "_" + m + "_meta_het_p")

        out.write("\n")

        next_var = get_next_variant(studs)
        if not args.quiet:
            print("NEXT VARIANTS")
            for v in next_var:
                print(v)
        matching_studies = [(studs[i],v) for i,v in enumerate(next_var) if v is not None]

        while len(matching_studies)>0:

            d = matching_studies[0][1]
            outdat = [ d.chr, d.pos, d.ref, d.alt]
            v = "{}:{}:{}:{}".format(*outdat)
            outdat.append(v)

            for i,_ in enumerate(studs):
                if next_var[i] is not None:
                    outdat.extend([format_num(next_var[i].beta), format_num(next_var[i].se), format_num(next_var[i].pval) ])
                    outdat.extend([ c for c in next_var[i].extra_cols ])

                    # meta analyse pairwise only with the leftmost study
                    if not args.pairwise_with_first or i==0:
                        continue

                    if next_var[0] is not None:
                        met = do_meta( [(studs[0],next_var[0]), (studs[i],next_var[i])], methods=methods, is_het_test=False)
                        for m in met:
                            outdat.extend(m)
                    else:
                        outdat.extend(["NA"] * len(methods) * 4)
                else:
                    outdat.extend(['NA']  * (3 + len(studs[i].extra_cols) + (len(methods)*4 if args.pairwise_with_first and i>0 else 0) ) )

            outdat.append( str(len(matching_studies)) )

            met = do_meta( matching_studies, methods=methods, is_het_test=args.is_het_test )
            for m in met:
                if m is not None:
                    outdat.extend(m)
                else:
                    outdat.extend(['NA'] * n_meta_cols)

            if args.leave_one_out:
                for s,_ in enumerate(studs):
                    matching_studies_loo = [(studs[i], var) for i,var in enumerate(next_var) if s != i and var is not None]
                    outdat.append( str(len(matching_studies_loo)) )
                    if len(matching_studies_loo) > 0:
                        met = do_meta( matching_studies_loo, methods=methods, is_het_test=args.is_het_test )
                        for m in met:
                            if m is not None:
                                outdat.extend(m)
                            else:
                                outdat.extend(['NA'] * n_meta_cols)
                    else:
                        outdat.extend(['NA'] * n_meta_cols * len(methods))

            out.write( "\t".join([ str(o) for o in outdat]) + "\n" )

            next_var = get_next_variant(studs)
            if not args.quiet:
                print("NEXT VARIANTS")
                for v in next_var:
                    print(v)
            matching_studies = [(studs[i],v) for i,v in enumerate(next_var) if v is not None]

    subprocess.run(["bgzip","--force",args.path_to_res])
    subprocess.run(["tabix","-s 1","-b 2","-e 2",args.path_to_res + ".gz"])


if __name__ == '__main__':
    run()
