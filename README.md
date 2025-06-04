# Bioinformartics-Training
Implementation of common algorithms for sequence analysis - quality control, assembly, alignment

1) Base Qualities.py - A Python script for basic quality control and sequence composition analysis of FASTQ files. 

- Reads FASTQ files and extracts sequences and Phred+33 quality scores

- Generates a quality score histogram to visualize overall read quality

- Computes GC content by position across all reads and plots it

- Counts base frequencies (A, T, G, C, N) across all sequences


2) Assembly.py -  A Python implementation of greedy shortest common superstring (SCS) genome assembly using read overlaps. 

- Reads FASTQ files and extracts sequences

- Finds maximal pairwise overlaps between reads using a k-mer index

- Greedily merges reads to assemble a genome-like superstring

- Reports assembly statistics, such as total length and base counts

This is a simplified assembler that could be useful for educational purposes / small-scale synthetic datasets.

3) Genome Search - Naive Matching.py - A toolkit for naive read alignment and search in DNA sequences.

This script provides a foundational framework for understanding how read mapping and pattern matching work in bioinformatics, using:

- Exact and approximate pattern matching

- Support for reverse complements

- Synthetic and real read input

- Match reporting and summary statistics

##### Example Use Cases
- Generate artificial reads from a viral genome and align them exactly

- Load real Illumina reads and check for exact or reverse complement matches

- Search for short patterns allowing mismatches (e.g. mutation-tolerant motifs)

- This script is ideal for learning and testing basic DNA sequence alignment methods without relying on complex tools like BWA or Bowtie.



