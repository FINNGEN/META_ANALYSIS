#! /usr/bin/env python3

import argparse
import gzip
from functools import partial
import sys
import subprocess
import uuid
import os

parser = argparse.ArgumentParser(description='Add 38 positions to summary file')
parser.add_argument("file", help=" Whitespace separated file with either single column giving variant ID in chr:pos:ref:alt or those columns separately")

parser.add_argument("-var", help="Variant column in chr:pos:ref:alt.")
parser.add_argument("-chr", help="Chromosome column name")
parser.add_argument("-pos", help="Position column name")
parser.add_argument("-ref", help="ref column name")
parser.add_argument("-alt", help="alt column name")
parser.add_argument("-chain_file", help="optional chain file path")
parser.add_argument("-no_clean", action="store_true" , default=False, help="Do not clean intermediate files after execution")


args = parser.parse_args()

script_dir =os.path.dirname(os.path.abspath(__file__))

if not args.chain_file:
    CHAINFILE= script_dir+ "/../data/hg19ToHg38.over.chain.gz"
else:
    CHAINFILE=args.chain_file

print("CHAINFILE" + CHAINFILE)
liftOver="liftOver"
chr = None
pos = None
ref = None
alt = None

get_dat_func = None
joinsortargs = ["-v variant"]

def get_dat_var(line, index):
    d = line[index].split(":")
    if len(d)<4:
        print("WARNING: Not properly formatted variant id in line: " + line, file=sys.stderr )
        return None
    return d

openf = open

suffix = args.file.split(".")[-1]

if suffix.lower() in ["gz","zip","bgzip"]:
    openf= gzip.open

with openf( args.file ,'rt') as res:
    header = res.readline().rstrip("\n").split()

    if args.var is None:
        if args.chr is None or args.pos is None or args.ref is None or args.alt is None:
            raise Exception("If var column not specified you must specify -chr -pos -ref and -alt")

        cols = [args.chr, args.pos, args.ref, args.alt]
        exists = [v in header for v in cols ]

        def f(d):
            i,v = d
            return not v

        if not all(exists):
            raise Exception("All required column names not in given data. Missing columns:" + " ".join([ c for c in [ cols[i] for (i,v) in filter(f, enumerate(exists)) ] ]) )

        chr = header.index(args.chr)
        pos = header.index(args.pos)
        ref = header.index(args.ref)
        alt = header.index(args.alt)
        joinsortargs = ["--chr",chr+1,"--pos",pos+1,"--ref", ref+1, "--alt", alt+1]
        get_dat_func = lambda line:  (line[chr], line[pos], line[ref], line[alt])
    else:
        var = header.index(args.var)
        joinsortargs = ["--var",var+1]
        get_dat_func = partial(get_dat_var,index=var)
    tempbed = str(uuid.uuid4())

    f = os.path.basename(args.file)
    tmpbed = f + tempbed + "_tmp.bed"
    with open(tmpbed, 'w') as bed:
        for line in res:
            vardat = get_dat_func(line.rstrip("\n").split())
            try:
                bed.write( "{}\t{}\t{}\t{}".format(vardat[0] if vardat[0].startswith('chr') else "chr"+vardat[0], str(int(vardat[1])-1), str(int(vardat[1]) + max(len(vardat[2]),len(vardat[3]) ) -1), ":".join([vardat[0],vardat[1],vardat[2],vardat[3]])) + "\n" )
            except ValueError:
                print("Ignoring unexpected chromosome position: " + ":".join([vardat[0],vardat[1],vardat[2],vardat[3]]))
                pass

    temp = f + str(uuid.uuid4())
    subprocess.run([liftOver, tmpbed ,CHAINFILE, temp + "_lifted", temp + "_errors"])

    joincmd = [script_dir +"/joinsort.sh", args.file, temp + "_lifted" ]
    joincmd.extend([ str(v) for v in joinsortargs])
    subprocess.run(joincmd)

    if not args.no_clean:
        print("Deleting intermediate files " + temp + "_lifted and " +tmpbed)
        subprocess.run(["rm",temp + "_lifted"])
        subprocess.run(["rm",tmpbed])
