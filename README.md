# Greedy_shortest_common_superstring
Context: By using and modifying Python programs provided in the Coursera course entitled 'Algorithms for DNA Sequencing', I used the shortest common superstring algorithm, the greedy shortest common superstring algorithm, and an accelerated version of the greedy shortest common superstring algorithm to assemble sequencing reads.

About:
1. overlap(a, b, min_length) finds the length of the longest suffix of a which overlaps with a prefix of b.
2. scs(ss) implements the shortest common superstring algorithm to find the set of shortest common superstrings of given strings ss.
3. Given a set of reads, pick_maximal_overlap finds the pair which overlap the most and calculates the length of the overlap.
4. greedy_scs implements the greedy shortest common superstring algorithm to find the set of shortest common superstrings of given strings.
5. pick_maximal_overlap_index and greedy_scs_index are accelerated versions of 3 and 4 respectively. Acceleration is achieved by building an k-mer index so that not every permutation of reads is considered.
6. Sequencing errors, ploidy, or reads coming from the reverse strand are not considered by the program. Also, it does not consider the effects of repeats.

Files:
1. main.py and ads1_week4_reads.fq must be in the same directory.
2. main.py should be implemented in Python 3.7.
3. ads1_week4_reads.fq is a FASTQ file containing synthetic sequencing reads from a mystery virus.

Modules:
1. itertools
