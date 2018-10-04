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

flip = {"A":"T","C":"G","T":"A","G":"C"}

def flip_strand( allele):
    return "".join([ flip[a] for a in allele])


class MetaDat:

    def __init__(self, chr, pos, ref, alt, beta, pval, se=None):
        self.chr = chr
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.beta = beta
        self.pval = pval
        self.se = float(se) if se is not None else None

    def __eq__(self, other):

        return self.chr == other.chr and self.pos == other.pos and self.ref == other.ref and self.alt == other.alt

    def equalize(self, other:'MetaDat') -> bool:
        """
            Checks if this metadata is the same variant (possibly different strand or ordering of alleles)
            returns: true if the same (flips effect direction and ref/alt alleles if necessary) or false if not the same variant
        """

        if (self.chr == other.chr and self.pos == other.pos):
            flip_ref =  flip_strand(other.ref)
            flip_alt =  flip_strand(other.alt)

            if( (self.ref == other.ref or self.ref==flip_ref) and (self.alt == other.alt or self.alt == flip_alt)):
                return True
            elif (self.ref == other.alt or self.ref == flip_alt) and (self.alt == other.ref or self.alt==flip_alt ) :
                self.beta = -1 * self.beta
                t = self.alt
                self.alt = self.ref
                self.ref = t
                return True
            else:
                return False
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

    def __init__(self, conf):
        self.conf =conf
        self.future = None

        for v in Study.REQUIRED_CONF:
            if v not in self.conf:
                raise Exception("Meta configuration for study must contain required elements: "
                    + ",".join(Study.REQUIRED_CONF.keys() ) + ". Offending configuration: " + str(s))

            try:
                self.conf[v] = Study.REQUIRED_CONF[v](self.conf[v])
            except Exception as e:
                raise Exception("Illegal data type in configuration for field " + s[v] +
                    " in configuration: " + str(s) + ". ERR:" + str(e))

        for v in Study.OPTIONAL_FIELDS:
            if v not in self.conf:
                continue
            try:
                self.conf[v] = Study.OPTIONAL_FIELDS[v](self.conf[v])
            except Exception as e:
                raise Exception("Illegal data type in configuration for field " + s[v] +
                    " in configuration: " + str(s) + ". ERR:" + str(e))

        self.conf["fpoint"] = gzip.open(conf["file"],'rt')
        header = conf["fpoint"].readline().rstrip().split("\t")
        self.conf["h_idx"] = { k:header.index( self.conf[k] ) for k in Study.REQUIRED_DATA_FIELDS.keys() }

        for f in Study.OPTIONAL_FIELDS.keys():
            if f in self.conf:
                 self.conf["h_idx"][f] = header.index(self.conf[f])

    @property
    def n_cases(self):
        return self.conf["n_cases"]

    @property
    def n_controls(self):
        return self.conf["n_cases"]

    @property
    def name(self):
        return self.conf["name"]


    def has_std_err(self):
        return "se" in self.conf

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

        if not chr.endswith("_alt") and not chr.endswith("_random"):
            chr = chrord[chr]

        return MetaDat(chr,int(pos),ref,alt, eff, pval, se)

    def get_match(self, dat: MetaDat) -> MetaDat:
        """
            Reads current study until variant in 'dat' is reached or overtaken in chr pos orded.
            IF matching variant found (can flip alleles) the matching MetaDat(effect flipped if alleles flipped) is returned.
            input:
                dat: the variant to look for
            output: matching MetaDat in this study or None if no match.
        """

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




def get_studies(conf:str) -> List[Study]:
    """
        Reads json configuration and returns studies in the meta
    """

    studies_conf = json.load(open(conf,'r'))
    std_list = studies_conf["meta"]

    return [ Study(s) for s in studies_conf["meta"]]

def do_meta(study_list: List[ Tuple[Study, MetaDat]] ) -> Tuple[float, float] :
    '''
        Computes meta-analysis between all studies and data given in the std_list
        input:
            study_list: studies and data in tuples
        output:
            tuple in which first element is effective sample size weighted meta and second is std err weighted meta. Second element is none if
            all studies did not have optional std err defined
    '''

    effs_size = []
    tot_size = 0

    effs_se = []
    tot_se = 0
    for s in study_list:
        study = s[0]
        dat = s[1]
        eff_size = ( (4 * study.n_cases *  study.n_controls  ) / ( study.n_cases+  study.n_controls ))
        zscore = math.copysign(1, dat.beta) * math.sqrt(chi2.isf(dat.pval, df=1))
        
        if dat.se is not None:
            effs_se.append( (1/ dat.se )* zscore )
            tot_se+= (1/ (dat.se * dat.se) )

        effs_size.append( math.sqrt(eff_size) * zscore)
        tot_size+=eff_size

    size_meta_p = scipy.stats.norm.sf( abs( sum( effs_size ) ) / math.sqrt(tot_size) )
    stderr_meta_p = scipy.stats.norm.sf( abs( sum( effs_se ) ) / math.sqrt(tot_se) ) if len(effs_se)==len(study_list) else None
    return (size_meta_p,stderr_meta_p)


def run():
    '''
        This module generates matrix from external single association results for fast access to browsingself.
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

    parser = argparse.ArgumentParser(description="Create tabixed big matrix for external results")
    parser.add_argument('config_file', action='store', type=str, help='Configuration file ')
    parser.add_argument('path_to_res', action='store', type=str, help='Result file')

    args = parser.parse_args()

    studs = get_studies(args.config_file)

    outfile = args.path_to_res

    with open( outfile, 'w' ) as out:

        out.write("\t".join(["CHR","POS","REF","ALT", studs[0].name + "_beta", studs[0].name + "_pval"  ]))

        ## align to leftmost STUDY
        for oth in studs[1:len(studs)]:
            out.write( "\t" +  "\t".join( [ oth.name + "_beta", oth.name + "_pval", studs[0].name + "_" + oth.name +"_meta_p"] ))
            if studs[0].has_std_err() and oth.has_std_err():
                out.write("\t" + studs[0].name + "_" + oth.name +"se_meta_p")

        if sum( map( lambda x: x.has_std_err(), studs ))>1:
            out.write("\tall_meta_N\tall_meta_p\tall_se_meta_p\n")
        else:
            out.write("\tall_meta_N\tall_meta_p\n")


        while True:
            ## change to class based get_data stuff...
            d = studs[0].get_next_data()

            if d is None:
                break
            matching_studies = [ (studs[0],d) ]
            outdat = [ d.chr, d.pos, d.ref, d.alt, numpy.format_float_scientific(d.beta, precision=2), numpy.format_float_scientific(d.pval, precision=2)  ]

            for oth in studs[1:len(studs)]:
                match_dat = oth.get_match(d)
                if match_dat is not None:
                    matching_studies.append( (oth,match_dat) )
                    met = do_meta( [(studs[0],d), (oth,match_dat)] )
                    outdat.extend([numpy.format_float_scientific(match_dat.beta, precision=2), numpy.format_float_scientific(match_dat.pval, precision=2), numpy.format_float_scientific(met[0], precision=2) ])

                    if( studs[0].has_std_err() & oth.has_std_err() ):
                        outdat.append( numpy.format_float_scientific(met[0], precision=2)  )
                else:
                    outdat.extend(['NA']  * 3)

            if len(matching_studies)>1:
                met = do_meta( matching_studies )
                outdat.append( str(len(matching_studies)) )
                outdat.append( numpy.format_float_scientific(met[0], precision=2) )
                if all( map( lambda x: x.has_std_err(), studs )):
                    outdat.append( numpy.format_float_scientific(met[1], precision=2) )

            else:
                outdat.extend(["NA"] * (2 + (1 if sum( map( lambda x: x.has_std_err(), studs ))>1 else 0) ))

            out.write( "\t".join([ str(o) for o in outdat]) + "\n" )

if __name__ == '__main__':
    run()
