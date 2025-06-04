# Bioinformatics Training

Implementation of common algorithms for sequence analysis including quality control, assembly, and alignment.

---

## 1) Base Qualities.py

A Python script for basic quality control and sequence composition analysis of FASTQ files.

- Reads FASTQ files and extracts sequences and Phred+33 quality scores  
- Generates a quality score histogram to visualize overall read quality  
- Computes GC content by position across all reads and plots it  
- Counts base frequencies (A, T, G, C, N) across all sequences  

**Testing data:**  
Source FASTQ file `ERR037900_1.first1000.fastq` from Coursera Genomics Specialisation.

---

## 2) Assembly.py

A Python implementation of greedy shortest common superstring (SCS) genome assembly using read overlaps.

- Reads FASTQ files and extracts sequences  
- Finds maximal pairwise overlaps between reads using a k-mer index  
- Greedily merges reads to assemble a genome-like superstring  
- Reports assembly statistics such as total length and base counts  

**Testing data:**  
Custom synthetic short reads file `ads1_week4_reads.fq` sourced from a module in Coursera Bioinformatics Specialisation.

---

## 3) Genome Search - Naive Matching.py

A toolkit for naive read alignment and search in DNA sequences.

Provides a foundational framework for understanding read mapping and pattern matching using:

- Exact and approximate pattern matching  
- Support for reverse complements  
- Synthetic and real read input  
- Match reporting and summary statistics  

**Testing data:**  
Custom short patterns (`ERR266411_1.first1000.fastq`) matched to the Lambda virus genome reference (`lambda_virus.fa`).

**Example use cases:**  
- Generate artificial reads from a viral genome and align them exactly  
- Load real Illumina reads and check for exact or reverse complement matches  
- Search for short patterns allowing mismatches (e.g. mutation-tolerant motifs)  

This script can be useful for learning and testing basic DNA sequence alignment methods without relying on complex tools like BWA or Bowtie.

---

## 4) Motif-Search.py

Implements multiple motif-finding algorithms for bioinformatics applications:

- **Gibbs Sampling**  
- **Randomized Motif Search**  
- Utilities for scoring, profile creation, and k-mer selection  

**Testing data:**  
Custom artificially created DNA sequences.

---

## Usage (for all scripts)

```bash
python <name-of-the-file>.py
