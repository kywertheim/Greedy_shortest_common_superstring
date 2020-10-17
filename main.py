# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 16:45:10 2020

@author: Kenneth
"""

"""
Import modules needed.
"""
import itertools

"""
This function reads a FASTQ file.
"""
def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

"""
This function finds the length of the longest suffix of a which overlaps with a prefix of b.
"""
def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

"""
This function finds the set of shortest common superstrings of given strings.
Note that the given strings must have the same length.
"""
def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = []
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if len(shortest_sup) == 0 or len(sup) < len(shortest_sup[0]):
            shortest_sup = [sup]  # found shorter superstring
        elif len(sup) == len(shortest_sup[0]):
            shortest_sup.append(sup)
    return shortest_sup  # return shortest

"""
Given a set of reads, this function finds the pair which overlap the most and calculates the length of the overlap.
"""
def pick_maximal_overlap(reads, k):
    reada, readb = None, None
    best_olen = 0
    for a,b in itertools.permutations(reads, 2):
        olen = overlap(a, b, k)
        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen
    return reada, readb, best_olen

"""
This function implements the greedy shortest common superstring algorithm.
"""
def greedy_scs(reads, k):
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        print(len(reads))
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    return ''.join(reads)

"""
This is an accelerated version of pick_maximal_overlap(reads, k).
This is achieved by building an k-mer index so that not every permutation of reads is considered.
"""
def pick_maximal_overlap_index(reads, k):
    index = {}
    for read in reads:
        kmers = []
        for i in range(len(read) - k + 1):
            kmers.append(read[i:i+k])
        for kmer in kmers:
            if kmer not in index:
                index[kmer] = set()
            index[kmer].add(read)
    for read in reads:
        for i in range(len(read)-k+1):
            dummy = read[i:i+k]
            if dummy not in index:
                index[dummy] = set()
            index[dummy].add(read)
    reada, readb = None, None
    best_olen = 0
    for a in reads:
        for b in index[a[-k:]]:
            if a != b:
                olen = overlap(a, b, k)
                if olen > best_olen:
                    reada, readb = a, b
                    best_olen = olen
    return reada, readb, best_olen

"""
This function implements the greedy shortest common superstring algorithm using an accelerated version of pick_maximal_overlap(reads, k).
"""
def greedy_scs_index(reads, k):
    read_a, read_b, olen = pick_maximal_overlap_index(reads, k)
    while olen > 0:
        print(len(reads))
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap_index(reads, k)
    return ''.join(reads)      

"""
Load the FASTQ file containing synthetic sequencing reads from a mystery virus.
""" 
seq, qual = readFastq('ads1_week4_reads.fq')

"""
Question 1.
Assemble the given reads and report the length of the shortest common superstrings.
"""
ss = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
print(len(scs(ss)[0]))

"""
Question 2.
Assemble the given reads and report the number of shortest common superstrings.
"""
ss = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
print(len(scs(ss)))

"""
Question 3.
Assemble the synthetic sequencing reads to form the mystery virus genome.
How many As are there in the full, assembled genome?
greedy_scs_index is more efficient than greedy_scs.
"""
genome = greedy_scs(seq, 10)
genome = greedy_scs_index(seq, 10)
print(genome.count('A'))

"""
Question 4.
Assemble the synthetic sequencing reads to form the mystery virus genome.
How many Ts are there in the full, assembled genome?
greedy_scs_index is more efficient than greedy_scs.
"""
genome = greedy_scs(seq, 10)
genome = greedy_scs_index(seq, 10)
print(genome.count('T'))