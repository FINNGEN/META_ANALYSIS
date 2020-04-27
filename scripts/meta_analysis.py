#!/usr/bin/env python3
import argparse
import json
import gzip
from collections import namedtuple
import sys
import math
from scipy.stats import chi2
import scipy.stats
import numpy
from typing import Dict, Tuple, List
import subprocess
from collections import deque
import re

chrord = { "chr"+str(chr):int(chr) for chr in list(range(1,23))}
chrord["X"] = 23
chrord["Y"] = 24
chrord["MT"] = 25
chrord["chrX"] = 23
chrord["chrY"] = 24
chrord["chrMT"] = 25
chrord.update({str(chr):int(chr) for chr in list(range(1,25)) } )

re_allele = re.compile('^[ATCG]+$', re.IGNORECASE)


def het_test( effs_sizes, weights, effs_size_meta):
    k=len(effs_sizes)

    effs_sizes_array=numpy.array(effs_sizes)
    weights_array=numpy.array(weights)    
    eff_dev=weights_array*((effs_sizes_array-effs_size_meta)**2)
    sum_eff_dev=numpy.sum(eff_dev)

    return scipy.stats.distributions.chi2.sf(sum_eff_dev, k-1)

def n_meta( studies : List[Tuple['Study','VariantData']], is_het_test = False ):
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
    if is_het_test:
        het_p=het_test(effs_size_org, weights, beta_meta)
    else:
        het_p=None
    #TODO se
    return ( beta_meta, None, max(sys.float_info.min * sys.float_info.epsilon, 2 * scipy.stats.norm.sf( abs( sum( effs_size ) ) / math.sqrt(tot_size) )), het_p)


def inv_var_meta( studies : List[Tuple['Study','VariantData']], is_het_test = False):

    weights = []
    effs_size_org = []

    effs_inv_var = []
    sum_inv_var=0
    for s in studies:
        study = s[0]
        dat = s[1]
        if dat.se is None or dat.se==0:
            print("Standard error was none/zero for variant " + str(dat) + " in study " + study.name, file=sys.stderr)
            break
        var = (dat.se * dat.se)

        inv_var =  (1/var)
        sum_inv_var+=inv_var
        effs_inv_var.append( inv_var *  dat.beta )

        weights.append(inv_var)
        effs_size_org.append(dat.beta)

    beta_meta=sum(effs_inv_var)/ sum_inv_var
    if is_het_test:
        het_p=het_test(effs_size_org, weights, beta_meta)
    else:
        het_p=None
    return (beta_meta, math.sqrt(1/sum_inv_var), max(sys.float_info.min * sys.float_info.epsilon, 2 * scipy.stats.norm.sf(abs(sum(effs_inv_var) / math.sqrt(sum_inv_var) ))), het_p) if len(effs_inv_var)==len(studies) else None


def variance_weight_meta( studies : List[Tuple['Study','VariantData']], is_het_test = False ):

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
            break
        weight =  (1/dat.se) * dat.z_score
        sum_weights+=weight
        sum_betas+= weight * dat.beta
        effs_se.append( weight * numpy.sign(dat.beta)  )
        tot_se+=1/ (dat.se * dat.se)
        
        weights.append(weight)
        effs_size_org.append(dat.beta)

    beta_meta=sum_betas / sum_weights
    if is_het_test:
        het_p=het_test(effs_size_org, weights, beta_meta)
    else:
        het_p = None
    #TODO SE
    return (beta_meta, None, max(sys.float_info.min * sys.float_info.epsilon, 2 * scipy.stats.norm.sf( abs( sum( effs_se ) ) /  math.sqrt(tot_se))), het_p) if len(effs_se)==len(studies) else None
    

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

    def is_equal(self, other:'VariantData') -> bool:
        """
            Checks if this VariantData is the same variant (possibly different strand or ordering of alleles)
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

    def equalize_to(self, other:'VariantData') -> bool:
        """
            Checks if this VariantData is the same variant as given other variant (possibly different strand or ordering of alleles)
            If it is, changes this variant's alleles and beta accordingly
            returns: true if the same (flips effect direction and ref/alt alleles if necessary) or false if not the same variant
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
                    t = self.alt
                    self.alt = self.ref
                    self.ref = t
                    return True

            elif (self.ref == other.alt and self.alt == other.ref) :
                self.beta = -1 * self.beta if self.beta is not None else None
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
                self.ref =flip_strand(self.alt)
                self.alt = flip_strand(self.ref)
                return True

        return False

    @property
    def z_score(self):
        '''
            Lazy compute unsigned z-score
        '''
        if self.z_scr is None:
            self.z_scr = math.sqrt(chi2.isf(self.pval, df=1))
        return self.z_scr

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

    def __init__(self, conf, chrom=None, dont_allow_space=False):
        '''
        chrom: a chromosome to limit to or None if all chromosomes
        dont_allow_space: boolean, don't treat space as field delimiter (only tab)
        '''
        self.conf =conf
        self.chrom = chrom
        self.dont_allow_space = dont_allow_space
        self.future = deque()
        self.eff_size= None
        self.z_scr = None
        self.prev_var = None
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
        if self.dont_allow_space:
            header = conf["fpoint"].readline().rstrip().split('\t')
        else:
            header = conf["fpoint"].readline().rstrip().split()

        for k in Study.REQUIRED_DATA_FIELDS.keys():
            if self.conf[k] not in header:
                raise Exception("Required headers not in data in study " + self.conf["name"] + ". Missing:" + ",".join([ self.conf[k] for k in Study.REQUIRED_DATA_FIELDS.keys() if self.conf[k] not in header])  )
        self.conf["h_idx"] = { k:header.index( self.conf[k] ) for k in Study.REQUIRED_DATA_FIELDS.keys() }

        for f in Study.OPTIONAL_FIELDS.keys():
            if f in self.conf:
                 if self.conf[f] not in header:
                     raise Exception("Configured column " + self.conf[f] + " not found in the study results " + self.conf["name"])
                 self.conf["h_idx"][f] = header.index(self.conf[f])

        if "extra_cols" in self.conf:
            for c in self.conf["extra_cols"]:
                if c not in header:
                    raise Exception("Configured column " + self.conf[c] + " not found in the study results " + self.conf["name"])
                self.conf["h_idx"][c] = header.index(c)
        else:
             self.conf["extra_cols"] = []



    @property
    def n_cases(self):
        return self.conf["n_cases"]

    @property
    def n_controls(self):
        return self.conf["n_cases"]

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

    def get_next_data(self, just_one =False) -> List[VariantData]:
        """
            Returns a list of variants. List containts >1 elements if they are on the same position and just_one ==False.
            args:
                just_one: always returns only the next variant in order and not all next with the same position
            returns: list of next variants
        """
        if len(self.future)>0:
            ## only return variants with same position so that possible next variant position stored stays
            f = [ (i,v) for i,v in enumerate(self.future) if i==0 or (v.chr==self.future[i-1].chr and  v.pos==self.future[i-1].pos) ]
            for i,v in reversed(f):
                 del self.future[i]
            return [ v for i,v in f ]

        vars = list()
        while True:
            chr = None
            ref = None
            alt = None
            l = None
            ## loop ignoring  alternate contigs and non-ATCG alleles for now.
            while chr is None or chr not in chrord or (self.chrom is not None and chr != self.chrom) or re_allele.match(ref) is None or re_allele.match(alt) is None:
                l = self.conf["fpoint"].readline()
                if l=="":
                    return None

                if self.dont_allow_space:
                    l = l.rstrip().split('\t')
                else:
                    l = l.rstrip().split()
                chr = l[self.conf["h_idx"]["chr"]].replace("chr", "")
                ref = l[self.conf["h_idx"]["ref"]]
                alt = l[self.conf["h_idx"]["alt"]]

            pos = l[self.conf["h_idx"]["pos"]]
            eff = l[self.conf["h_idx"]["effect"]]
            pval = l[self.conf["h_idx"]["pval"]]

            pos = int(float(pos))

            se = l[self.conf["h_idx"]["se"]] if "se" in self.conf["h_idx"] else None

            effect_type = self.conf["effect_type"]
            try:
                pval = float(pval)
                eff = float(eff)
            except Exception as e:
                pval = None
                eff = None

            if( effect_type=="or" and eff):
                eff = math.log(eff)

            chr = chrord[chr]
            extracols = [ l[self.conf["h_idx"][c]] for c in self.conf["extra_cols"] ]

            v = VariantData(chr,pos,ref,alt, eff, pval, se, extracols)

            if self.prev_var is not None and v < self.prev_var:
                raise Exception("Disorder in study " + self.conf['name'] + " in file " + self.conf['file'] + ". Sort all summary statistic files by chromosome and then position and rerun.\nOffending line: " + "\t".join(l))
            self.prev_var = v
            if len(vars)==0 or ( vars[0].chr == v.chr and vars[0].pos == v.pos  ):
                added=False
                for v_ in vars:
                    if v.is_equal(v_):
                        print('ALREADY ADDED FOR STUDY ' + self.name + ': ' + str(v))
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
            if v.equalize_to(dat):
                del otherdats[i]
                self.put_back(otherdats)
                return v

        ## no match but stayed in the same pos. add variants back to future queue
        self.put_back(otherdats)
        return None


    def put_back(self, VariantData):
        for m in VariantData:
            ## the future in next position will be always kept last
            self.future.appendleft(m)


def get_studies(conf:str, chrom, dont_allow_space) -> List[Study]:
    """
        Reads json configuration and returns studies in the meta
    """

    studies_conf = json.load(open(conf,'r'))
    std_list = studies_conf["meta"]

    return [ Study(s, chrom, dont_allow_space) for s in studies_conf["meta"]]

def do_meta(study_list: List[ Tuple[Study, VariantData]], methods: List[str], is_het_test) -> List[float] :
    '''
        Computes meta-analysis between all studies and data given in the std_list
        input:
            study_list: studies and data in tuples
        output:
            list of  tuples (effect_size, p-value) for each method in the same order as methods were given
    '''
    return [ SUPPORTED_METHODS[m](study_list, is_het_test) for m in methods ]

def format_num(num, precision=2):
    return numpy.format_float_scientific(num, precision=precision) if num is not None else "NA"

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

        for v in dats[i]:
            if v.equalize_to(first):
                res.append(v)
            else:
                s.put_back([v])
        if len(res)<i+1:
            res.append(None)

    return res



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

    parser = argparse.ArgumentParser(description="Run x-way meta-analysis")
    parser.add_argument('config_file', action='store', type=str, help='Configuration file ')
    parser.add_argument('path_to_res', action='store', type=str, help='Result file')

    parser.add_argument('methods', action='store', type=str, help='List of meta-analysis methods to compute separated by commas.'
            + 'Allowed values [n,inv_var,variance]', default="inv_var")

    parser.add_argument('--not_quiet', action='store_false', dest='quiet', help='Print matching variants to stdout')
    parser.set_defaults(quiet=True)

    parser.add_argument('--leave_one_out', action='store_true', help='Do leave-one-out meta-analysis')
    parser.set_defaults(leave_one_out=False)

    parser.add_argument('--is_het_test', action='store_true', help='Do heterogeneity tests based on Cochrans Q and output het_p')
    parser.set_defaults(het_test=False)

    parser.add_argument('--pairwise_with_first', action='store_true', help='Do pairwise meta-analysis with the first given study')
    parser.add_argument('--dont_allow_space', action='store_true', help='Do not allow space as field delimiter')

    parser.add_argument('--chrom', action='store', type=str, help='Restrict to given chromosome')

    args = parser.parse_args()

    studs = get_studies(args.config_file, args.chrom, args.dont_allow_space)

    methods = []

    for m in args.methods.split(","):
        if m not in SUPPORTED_METHODS:
            raise Exception("Unsupported meta method" + m + " given. Supported values" + ",".join(SUPPORTED_METHODS))
        methods.append(m)

    if "inv_var" in methods or "variance" in methods:
        for s in studs:
            if not s.has_std_err():
                raise Exception("Variance based method requested but not all studies have se column specified.")

    outfile = args.path_to_res

    with open( outfile, 'w' ) as out:

        out.write("\t".join(["#CHR","POS","REF","ALT","SNP", studs[0].name + "_beta", studs[0].name + "_sebeta", studs[0].name + "_pval"  ]))

        out.write( ("\t" if len(studs[0].extra_cols) else "") + "\t".join( [studs[0].name + "_" + c for c in studs[0].extra_cols] ) )
        ## align to leftmost STUDY
        for oth in studs[1:len(studs)]:
            out.write( "\t" +  "\t".join( [ oth.name + "_beta", oth.name + "_sebeta", oth.name + "_pval"] ))
            out.write( ("\t" if len(oth.extra_cols) else "") + "\t".join( [oth.name + "_" + c for c in oth.extra_cols] ) )

            if args.pairwise_with_first:
                for m in methods:
                    out.write("\t" + studs[0].name + "_" + oth.name + "_" +  m + "_meta_beta\t" + studs[0].name + "_" + oth.name + "_" +  m + "_meta_sebeta\t" + studs[0].name + "_" + oth.name + "_" +  m + "_meta_p")

        out.write("\tall_meta_N")
        for m in methods:
            if args.is_het_test:
                out.write("\tall_"+m+"_meta_beta\tall_"+m+"_meta_sebeta\tall_"+  m +"_meta_p\tall_"+ m +"_het_p")
            else:
                out.write("\tall_"+m+"_meta_beta\tall_"+m+"_meta_sebeta\tall_"+  m +"_meta_p")

        if args.leave_one_out:
            for s in studs:
                out.write("\t" + "leave_" + s.name + "_N")
                for m in methods:
                    out.write( "\t" +  "\t".join( ["leave_" + s.name + "_" + m + "_meta_beta", "leave_" + s.name + "_" + m + "_meta_sebeta", "leave_" + s.name + "_" + m + "_meta_p", "leave_" + s.name + "_" + m + "_meta_het_p"] ))

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
                        met = do_meta( [(studs[0],next_var[0]), (studs[i],next_var[i])], methods=methods , is_het_test=False)
                        for m in met:
                            outdat.append(format_num(m[0]))
                            outdat.append(format_num(m[1]))
                            outdat.append(format_num(m[2]))
                    else:
                        outdat.extend(["NA"] * len(methods) * 3)
                else:
                    outdat.extend(['NA']  * (3 + len(studs[i].extra_cols) + (len(methods)*3 if args.pairwise_with_first and i>0 else 0) ) )

            meta_res = []
            if len( matching_studies )>1:
                met = do_meta( matching_studies, methods=methods, is_het_test=args.is_het_test )
                for m in met:
                    meta_res.extend([format_num(num) for num in m[0:4]])
            else:
                meta_res.extend( [format_num(matching_studies[0][1].beta), format_num(matching_studies[0][1].se) , format_num(matching_studies[0][1].pval), 'NA']  * len(methods) )

            outdat.append( str(len(matching_studies)) )
            outdat.extend(meta_res)
            
            if args.leave_one_out:
                for s,_ in enumerate(studs):
                    matching_studies_loo = [(studs[i], var) for i,var in enumerate(next_var) if s != i and var is not None]
                    outdat.append( str(len(matching_studies_loo)) )
                    if len(matching_studies_loo) > 1:
                        met = do_meta( matching_studies_loo, methods=methods, is_het_test=args.is_het_test )
                        for m in met:
                            outdat.extend([format_num(num) for num in m[0:4]])
                    elif len(matching_studies_loo) == 1:
                        outdat.extend( [format_num(matching_studies_loo[0][1].beta), format_num(matching_studies_loo[0][1].se) , format_num(matching_studies_loo[0][1].pval), 'NA']  * len(methods) )
                    else:
                        outdat.extend(['NA'] * 4 * len(methods))

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
