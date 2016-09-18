#!/usr/bin/env python

from __future__ import print_function
import sys
from functools import partial
import string
import argparse
import pyximport
import numpy as np
pyximport.install(setup_args={'include_dirs': np.get_include()})
from parsing_pileup import (parseBases, qualToInt, processLine)


def getopt():
    parser = argparse.ArgumentParser(description='This program tries to convert samtools mipleup format to base count table with coverage threshold and quality threshold')
    parser.add_argument('-i','--input',help='mpileup format input filename (default: <->)',required=True)
    parser.add_argument('-c','--cov',default=10, type=int, help='coverage threshold (default: 10)')
    parser.add_argument('-q','--qual',default=33, type=int, help='quality threshold (default: 33)')
    args = parser.parse_args()
    return args

def main():
    # get arguments and call cython function
    args = getopt()
    qual_threshold = args.qual
    cov_threshold = args.cov
    filename = args.input
    header= ['chrom','start','end','ref_base','coverage','strand','A','C','T','G','deletions','insertions']
    print('\t'.join(header), file=sys.stdout)
    lineFunc = partial(processLine, qual_threshold, cov_threshold)
    handle = sys.stdin if filename == '-' else open(filename,'r')
    for lineno, line in enumerate(handle):
        lineFunc(line)
        if lineno % 100000 == 0:
            print('Parsed %i lines' %(lineno), file=sys.stderr)

if __name__ == '__main__':
    main()
