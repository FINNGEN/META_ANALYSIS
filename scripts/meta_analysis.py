#!/usr/bin/env python3
import argparse
import json
import gzip
from collections import namedtuple
import sys
import math
from scipy.stats import chi2
import scipy.stats

chrord = { "chr"+str(chr):int(chr) for chr in list(range(1,23))}
chrord["chrX"] = 23
chrord["chrT"] = 24
chrord["chrMT"] = 25
chrord.update({str(chr):int(chr) for chr in list(range(1,23)) } )

def check_eff_field(field):
    if field.lower() in ["beta","or"]:
        return field.lower()
    else:
        raise Exception("effect_type must be beta or OR")

required_conf = {"name":str,"file":str, "n_cases": int, "n_controls":int,
    "chr":str,"pos":str,"ref":str,"alt":str, "effect":str,
    "effect_type":check_eff_field,
    "pval":str}

required_in_meta_dat = {"chr":str,"pos":str,"ref":str,"alt":str, "effect":str,
    "pval":str}


flip = {"A":"T","C":"G","T":"A","G":"C"}

def flip_strand( allele):
    return [ flip[a] for a in allele]


class Study:

    def __init__(self, conf):
        self.conf =conf
        self.future = None

        self.REQUIRED_DATA_FIELDS = {"chr":str,"pos":str,"ref":str,"alt":str, "effect":str,
            "pval":str}

        self.REQUIRED_CONF = {"name":str,"file":str, "n_cases": int, "n_controls":int,
        "chr":str,"pos":str,"ref":str,"alt":str, "effect":str,
        "effect_type":check_eff_field,
        "pval":str}

        for v in self.REQUIRED_CONF:
            if v not in self.conf:
                raise Exception("Meta configuration for study must contain required elements: "
                    + ",".join(REQUIRED_CONF.keys() ) + ". Offending configuration: " + str(s))

            try:
                self.conf[v] = self.REQUIRED_CONF[v](self.conf[v])
            except Exception as e:
                raise Exception("Illegal data type in configuration for field " + s[v] +
                    " in configuration: " + str(s) + ". ERR:" + str(e))

        self.conf["fpoint"] = gzip.open(conf["file"],'rt')
        header = conf["fpoint"].readline().rstrip().split("\t")
        self.conf["h_idx"] = { k:header.index( self.conf[k] ) for k in self.REQUIRED_DATA_FIELDS.keys() }

    @property
    def n_cases(self):
        return self.conf["n_cases"]

    @property
    def n_controls(self):
        return self.conf["n_cases"]

    @property
    def name(self):
        return self.conf["name"]

    def get_next_data(self):

        if self.future is not None:
            f = self.future
            self.future =None
            return f

        l = self.conf["fpoint"].readline()
        if l=="":
            return None
        l = l.rstrip().split("\t")
        chr = l[self.conf["h_idx"]["chr"]]
        pos = l[self.conf["h_idx"]["pos"]]
        ref = l[self.conf["h_idx"]["ref"]]
        alt = l[self.conf["h_idx"]["alt"]]
        eff = l[self.conf["h_idx"]["effect"]]
        pval = l[self.conf["h_idx"]["pval"]]
        effect_type = self.conf["effect_type"]

        try:
            pval = float(pval)
            eff = float(eff)
        except Exception as e:
            pval = None
            eff = None

        if( effect_type=="or" and eff):
            eff = math.log(eff)

        if not chr.endswith("_alt") and not chr.endswith("_random"):
            chr = chrord[chr]

        return MetaDat(chr,int(pos),ref,alt, eff, pval)

    def get_match(self, dat):

        otherdat = self.get_next_data( )

        while otherdat is not None and otherdat.chr==dat.chr and otherdat.pos<dat.pos:
            otherdat = self.get_next_data()

        if otherdat is None:
            return None

        if otherdat.chr != dat.chr or otherdat.pos> dat.pos:
            self.put_back(otherdat)
            return None
        elif otherdat.equalize(dat):
            return otherdat

    def put_back(self, metadat):
        self.future =metadat

class MetaDat():

    def __init__(self, chr, pos, ref, alt, beta, pval):
        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.beta = beta
        self.pval = pval


    ('MetaDat','chr pos ref alt beta pval')
    def __eq__(self, other):

        return self.chr == other.chr and self.pos == other.pos and self.ref == other.ref and self.alt == other.alt

    def equalize(self, other):
        """
            Checks if this metadata is the same variant (possibly different strand or ordering of alleles)
            returns: true if the same (flips effect direction if necessary) or false if not the same variant
        """

        if (self.chr == other.chr and self.pos == other.pos):
            flip_ref =  flip_strand(other.ref)
            flip_alt =  flip_strand(other.alt)
            if( (self.ref == other.ref or self.ref==flip_ref) and (self.alt == other.alt or self.alt == flip_alt)):
                return True
            elif (self.ref == other.alt or self.ref == flip_alt) and (self.alt == other.ref or self.alt==flip_alt ) :
                self.beta = -1 + self.beta
                return True
            else:
                return False
    def __str__(self):
        return "chr:{} pos:{} ref:{} alt:{} beta:{} pval:{}".format(self.chr, self.pos, self.ref, self.alt, self.beta, self.pval)

def get_studies(conf):
    studies_conf = json.load(open(conf,'r'))
    std_list = studies_conf["meta"]

    return [ Study(s) for s in studies_conf["meta"]]

def do_meta(study_list):

    effs = []
    tot_size = 0

    ## check how to weight by MAF differences!!!
    ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5287121/

    for s in study_list:
        study = s[0]
        dat = s[1]
        eff_size = ( (4 * study.n_cases *  study.n_controls  ) / ( study.n_cases+  study.n_controls ))

        zscore = math.copysign(1, dat.beta) * math.sqrt(chi2.ppf(1-dat.pval, df=1))
        effs.append( math.sqrt(eff_size) * zscore)
        tot_size+=eff_size

    return 1-scipy.stats.norm.cdf( abs( sum( effs ) ) / math.sqrt(tot_size) )


def run():
    '''
        This module generates matrix from external single association results for fast access to browsingself.
        First parameter should be a path to a json configuration file with these elements:
            "name":"FINNGEN",
            "file":"/Users/mitja/projects/finngen/META_ANALYSIS/I9_AF.gz",
            "n_controls": 6570 ,
            "n_cases": 48378,
            "chr":"CHR",
            "pos":"POS",
            "ref":"Allele1",
            "alt":"Allele2",
            "effect":"BETA",
            "effect_type":"beta",
            "pval":"p.value"
        Second parameter should be a path to (empty/not existing) directory where the data should be stored
    '''

    parser = argparse.ArgumentParser(description="Create tabixed big matrix for external results")
    parser.add_argument('config_file', action='store', type=str, help='Configuration file ')
    parser.add_argument('path_to_res', action='store', type=str, help='Prefix filepath where the results will be saved')
    parser.add_argument('--chr', default="chr", action='store', type=str, help='chr column name in result files')
    parser.add_argument('--pos', default="pos", action='store', type=str, help='pos column name in result files')
    parser.add_argument('--ref', default="ref", action='store', type=str, help='ref column name in result files')
    parser.add_argument('--alt', default="alt", action='store', type=str, help='alt column name in result files. This is assumed to be the effective allele')
    parser.add_argument('--beta', default="beta", action='store', type=str, help='beta column name in result files')
    parser.add_argument('--pval', default="pval", action='store', type=str, help='pvalue column name in result files')
    parser.add_argument('--other_fields', action='store', type=str, help='comma separated list of other column names in result files')

    args = parser.parse_args()

    studs = get_studies(args.config_file)

    outfile = args.path_to_res

    with open( outfile, 'w' ) as out:

        out.write("\t".join(["CHR","POS","REF","ALT", studs[0].name + "_beta", studs[0].name + "_pval"  ]))

        ## align to leftmost STUDY
        for oth in studs[1:len(studs)]:
            out.write( "\t" +  "\t".join( [ oth.name + "_beta", oth.name + "_pval", studs[0].name + "_" + oth.name +"_meta_p"] ))
        out.write("\tall_meta_N\tall_meta_p\n")


        while True:
            ## change to class based get_data stuff...
            d = studs[0].get_next_data()

            if d is None:
                break
            matching_studies = [ (studs[0],d) ]
            outdat = [ d.chr, d.pos, d.ref, d.alt, d.beta, d.pval  ]

            for oth in studs[1:len(studs)]:
                match_dat = oth.get_match(d)
                if match_dat is not None:
                    matching_studies.append( (oth,match_dat) )
                    met = do_meta( [(studs[0],d), (oth,match_dat)] )
                    outdat.extend([str(match_dat.beta), str(match_dat.pval), str(met) ])
                else:
                    outdat.extend(['NA']  * 3)

            if len(matching_studies)>1:
                met = do_meta( matching_studies )
                outdat.append( str(len(matching_studies)) )
                outdat.append( str(met) )
            else:
                outdat.extend(["NA"] * 2 )

            out.write( "\t".join([ str(o) for o in outdat]) + "\n" )

if __name__ == '__main__':
    run()
