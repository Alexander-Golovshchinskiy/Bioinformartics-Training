import itertools

def readfastq(filename):
    """
    Reads a FASTQ file and returns a list of sequences and their qualities.
    
    Args:
        filename (str): Path to FASTQ file.
        
    Returns:
        sequences (list of str): List of DNA sequences.
        quality (list of str): List of quality strings corresponding to sequences.
    """
    sequences = []
    quality = []
    with open(filename, 'r') as fh:
        while True:
            fh.readline()  # Header line
            seq = fh.readline().rstrip()  # Sequence line
            fh.readline()  # Plus line
            qua = fh.readline().rstrip()  # Quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            quality.append(qua)
    return sequences, quality

def overlap(a, b, min_length=3):
    """
    Return length of longest suffix of 'a' matching a prefix of 'b' 
    that is at least 'min_length' characters long.
    
    Args:
        a (str): String a.
        b (str): String b.
        min_length (int): Minimum overlap length.
        
    Returns:
        int: Length of longest overlap or 0 if none found.
    """
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def build_kmer_index(reads, k):
    """
    Build an index mapping k-mers to reads containing them.
    
    Args:
        reads (list of str): List of reads.
        k (int): k-mer length.
        
    Returns:
        dict: Dictionary {kmer: set of reads containing kmer}
    """
    index = {}
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            if kmer not in index:
                index[kmer] = set()
            index[kmer].add(read)
    return index

def pick_maximal_overlap_index(reads, k):
    """
    Find pair of reads with maximal overlap >= k using k-mer index.
    
    Args:
        reads (list of str): List of reads.
        k (int): Minimum overlap length.
        
    Returns:
        tuple: (read_a, read_b, best_overlap_length) or (None, None, 0) if none found.
    """
    index = build_kmer_index(reads, k)
    reada, readb = None, None
    best_olen = 0
    for a in reads:
        suffix = a[-k:]
        for b in index.get(suffix, []):
            if a != b:
                olen = overlap(a, b, min_length=k)
                if olen > best_olen:
                    reada, readb = a, b
                    best_olen = olen
    return reada, readb, best_olen

def greedy_scs(reads, k):
    """
    Greedy shortest common superstring assembly of reads using maximal overlaps.
    
    Args:
        reads (list of str): List of reads.
        k (int): Minimum overlap length.
        
    Returns:
        str: Assembled superstring genome.
    """
    reads = reads.copy()
    read_a, read_b, olen = pick_maximal_overlap_index(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        merged = read_a + read_b[olen:]
        reads.append(merged)
        read_a, read_b, olen = pick_maximal_overlap_index(reads, k)
    return ''.join(reads)

def main():
    # Load reads from FASTQ file
    filename = 'ads1_week4_reads.fq'
    reads, _ = readfastq(filename)

    print(f"Number of reads loaded: {len(reads)}")

    # Set minimum overlap length (example: 30)
    k = 30

    # Run greedy assembly
    assembled_genome = greedy_scs(reads, k)

    # Print summary statistics
    print("Assembly completed.")
    print(f"Length of assembled genome: {len(assembled_genome)}")
    print(f"Count of 'A' in assembled genome: {assembled_genome.count('A')}")
    print(f"Count of 'T' in assembled genome: {assembled_genome.count('T')}")

if __name__ == '__main__':
    main()
