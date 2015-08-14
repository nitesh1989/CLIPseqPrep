"""
   python3 preprocess_py3.py [options] <fastq file>
guess the encoding of a stream of qual lines.
"""
import sys
import argparse
import subprocess as sub
#import tarfile
import gzip
import os

RANGES = {
    'phred33': (33, 74), ## Sanger, Illumina-1.8
    'phred64': (59, 104), ## Solexa, Illumina-1.3, Illumina-1.5
    # 'Sang Phred+33': (33, 73), ## Sanger
    # 'Solx Phred+64': (59, 104), ## Solexa
    # 'I1.3 Phred+64': (64, 104), ## Illumina-1.3
    # 'I1.5 Phred+64': (67, 104), ## Illumina-1.5
    # 'I1.8 Phred+33': (33, 74) ## Illumina-1.8
}


def get_qual_range(qual_str):
    """
    >>> get_qual_range("DLXYXXRXWYYTPMLUUQWTXTRSXSWMDMTRNDNSMJFJFFRMV")
    (68, 89)
    """

    vals = [ord(c) for c in qual_str]
    return min(vals), max(vals)

def get_encodings_in_range(rmin, rmax, ranges=RANGES):
    valid_encodings = []
    for encoding, (emin, emax) in ranges.items():
        if rmin >= emin and rmax <= emax:
            valid_encodings.append(encoding)
    return valid_encodings

def guess_encoding(nlines,filename):
    print("Reading qualities from fastq file...", file=sys.stderr)
    gmin, gmax  = 99, 0
    fastqfile = open(filename,'r')
    valid = []
    for i, line in enumerate(fastqfile.read(),start=1): # i=line number; line=sequence
        if i % 4 == 0:
            lmin, lmax = get_qual_range(line.rstrip()) # rstrip: default strips whitespace
            if lmin < gmin or lmax > gmax:
                gmin, gmax = min(lmin, gmin), max(lmax, gmax)
                valid = get_encodings_in_range(gmin, gmax)
                if len(valid) == 0:
                    print("no encodings for range: %s" % str((gmin, gmax)) + "\n",file=sys.stderr)
                    sys.exit()
                if len(valid) == 1 and nlines == -1:
                    return("\n".join(valid))
            if nlines > 0 and (i/4) > nlines:
                return("\n".join(valid)) #+"\n"+str((gmin, gmax)))
    return("\n".join(valid))

def call_prinseq(filename,phr):
    call = str.join("",("prinseq-lite -fastq ",filename," -params prinseq_params_",
        phr,".txt -out_good ",filename[:-6],"_good -out_bad ",filename[:-6],"_bad"))
    rc = os.system(call)
    return(rc)

def call_bowtie():
    return(0)

def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("-fastq", dest="file", help="Fastq formatted input file. By default,"
                 " reads from STDIN.", type=str)
    p.add_argument("-n","--nlines", dest="n", help="Number of qualily lines to test. By default,"
                 " tests until end of file or until it it possible to determine a single file-type (default: -1).",
                 type=int, default=-1) # dashes=optional args; no dash=required args
    args = p.parse_args()
    LOG = open("log.txt",mode="w")
    #fastqfile = gzip.open(args.file,mode='rt') ## 'rt' = read mode, text file
    phred = guess_encoding(args.n,args.file) ## get phred encoding value
    LOG.write("Successfully decoded fastq file quality scores as "+phred+".")

    retcode = call_prinseq(args.file,phred)
    if retcode == 0: LOG.write("\nReads successfully trimmed and filtered.")
    #retcode=subprocess.call(["prinseq-lite","-fastq",args.file,"-params",str.join("",("prinseq_params_",phred,".txt")),"-out_good","test_good","-out_bad","test_bad"])
    #info=subprocess.check_output(["gunzip","-c",args.file])
    #retcode=subprocess.call([call,"|","../../Tools/prinseq-lite-0.20.4/prinseq-lite","-fastq","stdin","-params",paramfile])

    retcode = call_bowtie()
    LOG.close()

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
