from __future__ import print_function
import string
import sys
from collections import Counter
from functools import partial
import numpy as np
cimport numpy as np
import re
complement_base = string.maketrans('actgnACTGN','TGACNTGACN')

def strandedBase(str bases_string):
    cdef:
        np.ndarray bases_negative
        np.ndarray bases_positive

    # extract lower and upper strand reads
    bases_negative = np.array(re.findall('[actg]',bases_string))
    bases_positive = np.array(re.findall('[ACTG]',bases_string))
    return bases_positive, bases_negative

cpdef int qualToInt(str q):
    cdef int qual = 0
    qual = ord(q) -33
    return qual

cpdef str qualityBases(str bases, str quals, int qual_threshold):
    cdef:
        np.ndarray bases_list
        np.ndarray quals_list
        str out_bases

    # extract high quality bases by numpy array
    bases_list = np.array(list(bases))
    quals_list = np.array(map(qualToInt,quals.strip()))
    assert len(bases_list)==len(quals_list), 'bases != quals ' +'\n'+ ''.join(bases) +'!\n' + ''.join(quals)
    quality_bases = bases_list[quals_list >= qual_threshold]
    out_bases = ''.join(quality_bases)
    return out_bases

cpdef int printLine(str chrom, int cov_threshold, str pos, str ref,
                    np.ndarray bases, str strand,
                    np.ndarray deletions, np.ndarray insertions):
    cdef:
        str bed_line_str
        int cov

    # print output line as bed format with base count and indel count
    bases_count = Counter(''.join(bases).upper())
    cov = np.sum(bases_count.values())
    bed_line = map(str, [chrom, pos, int(pos)+1, ref, cov, strand,
                         bases_count['A'], bases_count['C'], bases_count['T'], bases_count['G'],
                         len(deletions),len(insertions)])
    bed_line_str = '\t'.join(bed_line)
    print(bed_line_str,file=sys.stdout)
    return 0

def parseBases(str bases, str ref):
    cdef:
        int i = 0
        int j = 0
        int i_max = len(bases)
        int insertion = 0
        int deletion = 0
        str _bases = ''
        str insertion_bases = ''
        str deletion_bases = ''
        str c = ''

    # loop over the bases field
    # if it is ./, add to reference
    # skip for ^ or $ (start and end)
    # insertion/deletion if +/-, follwing number and bases are the indel bases
    # others are mapped base
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
            insert_digit = ''
            while bases[j].isdigit():
                insert_digit += bases[j]
                j += 1
            insert_count = int(insert_digit)
            i = j + insert_count - 1
            insertion_bases += bases[(j):i+1]
            insertion += insert_count
        elif c == '-':
            j = i + 1
            del_count = 0
            del_digit = ''
            while bases[j].isdigit():
                del_digit += bases[j]
                j += 1
            del_count = int(del_digit)
            i = j + del_count - 1
            deletion_bases += bases[(j):i+1]
            deletion += del_count
        elif c != '$':
            _bases += c
        i += 1
    return _bases, insertion, deletion, insertion_bases, deletion_bases

cpdef int processLine(int qual_threshold, int cov_threshold, str line):
    # define variables
    cdef:
        # define split input from mpileup line
        str chrom, pos, ref,cov, inbases, quals
        # define processed line result
        str bases, insertion_bases, deletion_bases
        int coverage, insertion, deletion
        np.ndarray bases_positive, bases_negative

    #split mpileup line
    fields = line.split('\t')
    chrom,  pos, ref, cov, inbases, quals = fields
    coverage = int(cov)
    quals = quals.strip()


    if coverage > cov_threshold:
        # using bases field to get information
        bases, insertion, deletion, insertion_bases, deletion_bases = parseBases(inbases, ref)
        assert len(bases) == len(quals),'Wrongly parsed!! ' + bases + ' ' + str(coverage)

        #extract high quality bases only
        bases = qualityBases(bases, quals, qual_threshold)

        # identify upper and lower strand reads
        bases_positive, bases_negative = strandedBase(bases)

        #identify indel from upper and lower strand
        insertion_positive, insertion_negative = strandedBase(insertion_bases)
        deletion_positive, deletion_negative = strandedBase(deletion_bases)

        # print line
        printFunc = partial(printLine, chrom, cov_threshold, pos)
        printed = map(printFunc, [ref.translate(complement_base), ref],
                                [bases_negative, bases_positive],
                                ['-','+'],
                                [deletion_negative, deletion_positive],
                                [insertion_negative, insertion_positive])
    return 0
