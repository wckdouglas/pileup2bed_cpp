#!/usr/bin/env python

from __future__ import print_function
import numpy
import sys
import re
from collections import Counter
import numpy as np
from functools import partial
import string
import argparse

complement_base = string.maketrans('actgnACTGN','TGACNTGACN')

def getopt():
    parser = argparse.ArgumentParser(description='This program tries to convert samtools mipleup format to base count table with coverage threshold and quality threshold')
    parser.add_argument('-i','--input',default='-',help='mpileup format input filename (default: <->)')
    parser.add_argument('-c','--cov',default=10, type=int, help='coverage threshold (default: 10)')
    parser.add_argument('-q','--qual',default=33, type=int, help='quality threshold (default: 33)')
    args = parser.parse_args()
    return args

def qualToInt(q):
    return ord(q) -33

def qualityBases(bases, quals, qual_threshold):
    bases = np.array(list(bases))
    quals_list = np.array(map(qualToInt,quals.strip()))
    #assert len(bases)==len(quals_list), 'bases != quals ' +'\n'+ ''.join(bases) +'!\n' + ''.join(quals)
    quality_bases = bases[quals_list >= qual_threshold]
    return ''.join(quality_bases)


def parseBases(bases, ref):
    i , i_max = 0, len(bases)
    _bases = ''
    insertion = 0
    deletion = 0
    insertion_bases, deletion_bases = '', ''
    while i < i_max:
        c = bases[i]
        if c == '.':
            _bases += ref
        elif c == ',':
            _bases += ref.translate(complement_base).lower()
        elif c == '^':
            i += 1
        elif c == '+':
            j = i + 1
            insert_count = 0
            while bases[j].isdigit():
                insert_count += int(bases[j])
                j += 1
            i = j + insert_count - 1
            insertion_bases += bases[(j):i+1]
            insertion += insert_count
        elif c == '-':
            j = i + 1
            del_count = 0
            while bases[j].isdigit():
                del_count += int(bases[j])
                j += 1
            i = j + del_count - 1
            deletion_bases += bases[(j):i+1]
            deletion += del_count
        elif c != '$':
            _bases += c
        i += 1
    return _bases, insertion, deletion, insertion_bases, deletion_bases

def printLine(chrom, cov_threshold, pos, ref, bases, strand, deletions, insertions):
    bases_count = Counter(''.join(bases).upper())
    cov = np.sum(bases_count.values())
    bed_line = map(str, [chrom, pos, int(pos)+1, ref, cov, strand,
                         bases_count['A'], bases_count['C'], bases_count['T'], bases_count['G'],
                         len(deletions),len(insertions)])
    bed_line = '\t'.join(bed_line)
    print(bed_line,file=sys.stdout)
    return 0

def strandedBase(bases_string):
    bases_negative = re.findall('[actg]',bases_string)
    bases_positive = re.findall('[ACTG]',bases_string)
    return bases_positive, bases_negative

def processLine(qual_threshold, cov_threshold, line):
    fields = line.split('\t')
    chrom,  pos, ref, cov, bases, quals = fields
    if int(cov) > cov_threshold:
        bases, insertion, deletion, insertion_bases, deletion_bases = parseBases(bases, ref)
        assert len(bases) == int(cov),'Wrongly parsed!! ' + bases + ' ' + cov
        bases = qualityBases(bases, quals, qual_threshold)
        bases_positive, bases_negative = strandedBase(bases)
        insertion_positive, insertion_negative = strandedBase(insertion_bases)
        deletion_positive, deletion_negative = strandedBase(deletion_bases)
        printFunc = partial(printLine, chrom, cov_threshold, pos)
        printed = map(printFunc, [ref.translate(complement_base), ref],
                                [bases_negative, bases_positive],
                                ['-','+'],
                                [deletion_negative, deletion_positive],
                                [insertion_negative, insertion_positive])
    return 0

def main():
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
        if lineno % 10000 == 0:
            print('Parsed %i lines' %(lineno), file=sys.stderr)

if __name__ == '__main__':
    main()
