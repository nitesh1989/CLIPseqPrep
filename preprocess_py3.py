#!/usr/bin/env python

"""
python3 preprocess_py3.py [options] <fastq file>
guess the encoding of a stream of qual lines.
"""
import sys
import argparse
import subprocess as sub
import gzip
import os
import ipdb
import pprint


RANGES = {
    'phred33': (33, 74),  # Sanger, Illumina-1.8
    'phred64': (59, 104),  # Solexa, Illumina-1.3, Illumina-1.5
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


def guess_encoding(nlines, filename):
    print("filename: {0}".format(filename))
    print("Reading qualities from fastq file...", file = sys.stderr)
    gmin, gmax = 99, 0
    if filename.endswith(".gz"):
        # The perl script is not built to handle fastq files with
        # .gz extension
        with gzip.open(filename, 'rt', encoding='utf-8') as f:
            fastqfile = f.read()
        new_filename = filename.replace(".gz", "")
        with open(new_filename, "w") as nf:
            nf.write(fastqfile)
    else:
        with open(filename, 'r', encoding='utf-8') as f:
            fastqfile = f.read()
#    ipdb.set_trace()
    valid = []
    for i, line in enumerate(fastqfile, start=1):  # i=line number; line=sequence
        if i % 4 == 0:
            lmin, lmax = get_qual_range(line.rstrip())  # rstrip: default strips whitespace
            if lmin < gmin or lmax > gmax:
                gmin, gmax = min(lmin, gmin), max(lmax, gmax)
                valid = get_encodings_in_range(gmin, gmax)
                if len(valid) == 0:
                    print("no encodings for range: %s" % str((gmin, gmax)) + "\n", file=sys.stderr)
                    sys.exit()
                if len(valid) == 1 and nlines == -1:
                    return("\n".join(valid))
            if nlines > 0 and (i / 4) > nlines:
                return("\n".join(valid))  # +"\n"+str((gmin, gmax)))
    return("\n".join(valid))


def call_prinseq(filename, phr):
    if filename.endswith(".gz"):
        filename = filename.replace(".gz", "")
    print("filename: {0}; phr: {1}".format(filename, phr))
    call = str.join("", ("prinseq-lite -fastq ", filename, " -params prinseq_params_", phr, ".txt -out_good ", filename[:-6], "_good -out_bad ", filename[:-6], "_bad"))
    rc = os.system(call)
    print("rc: %s" % rc)
    return(rc)


def call_bowtie():
    return(0)


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("-fastq", dest="file", help="Fastq formatted input file. By default,"
                 " reads from STDIN.", type=str)
    p.add_argument("-n", "--nlines", dest="n", help="Number of qualily lines to test. By default,"
                 " tests until end of file or until it it possible to determine a single file-type (default: -1).",
                 type=int, default=-1)  # dashes=optional args; no dash=required args
    args = p.parse_args()
    LOG = open("log.txt", mode="w")
    phred = guess_encoding(args.n, args.file)  # get phred encoding value
    LOG.write("Successfully decoded fastq file quality scores as " + phred + ".")

    retcode = call_prinseq(args.file, phred)
    if retcode == 0:
        LOG.write("\nReads successfully trimmed and filtered.")

    retcode = call_bowtie()
    LOG.close()

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()
