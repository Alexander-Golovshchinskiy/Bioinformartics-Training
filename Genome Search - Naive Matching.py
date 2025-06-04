import random

def read_genome(filename):
    """
    Read a genome sequence from a FASTA file (ignoring headers).
    
    Args:
        filename (str): Path to the FASTA file.
    
    Returns:
        genome (str): Genome sequence as a string.
    """
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome += line.strip()
    return genome

def naive_exact_search(pattern, text):
    """
    Naive exact matching of pattern in text.
    
    Args:
        pattern (str): Pattern string to search.
        text (str): Text string to search within.
    
    Returns:
        occurrences (list of int): Starting indices where pattern matches exactly.
    """
    occurrences = []
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i+len(pattern)] == pattern:
            occurrences.append(i)
    return occurrences

def naive_approx_search(pattern, text, max_mismatches):
    """
    Naive approximate matching allowing up to max_mismatches mismatches.
    
    Args:
        pattern (str): Pattern string to search.
        text (str): Text string to search within.
        max_mismatches (int): Maximum number of allowed mismatches.
    
    Returns:
        occurrences (list of int): Starting indices where pattern matches with <= max_mismatches.
    """
    occurrences = []
    plen = len(pattern)
    for i in range(len(text) - plen + 1):
        mismatches = 0
        for j in range(plen):
            if text[i+j] != pattern[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        if mismatches <= max_mismatches:
            occurrences.append(i)
    return occurrences

def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.
    
    Args:
        seq (str): DNA sequence.
    
    Returns:
        str: Reverse complement sequence.
    """
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join(complement[base] for base in reversed(seq))

def naive_search_with_reverse_complement(pattern, text):
    """
    Naive exact search allowing matches of pattern or its reverse complement.
    
    Args:
        pattern (str): Pattern string.
        text (str): Text string.
    
    Returns:
        occurrences (list of int): Unique starting indices matching pattern or its reverse complement.
    """
    comp_pattern = reverse_complement(pattern)
    occurrences = []
    plen = len(pattern)
    for i in range(len(text) - plen + 1):
        segment = text[i:i+plen]
        if segment == pattern or segment == comp_pattern:
            occurrences.append(i)
    return sorted(set(occurrences))

def generate_random_reads(genome, num_reads, read_len):
    """
    Generate random reads of given length from the genome.
    
    Args:
        genome (str): Reference genome string.
        num_reads (int): Number of reads to generate.
        read_len (int): Length of each read.
    
    Returns:
        reads (list of str): List of reads.
    """
    reads = []
    genome_len = len(genome)
    for _ in range(num_reads):
        start = random.randint(0, genome_len - read_len)
        reads.append(genome[start:start+read_len])
    return reads

def match_reads_to_genome(reads, genome, use_reverse_complement=False, max_mismatches=0):
    """
    Match a list of reads to the genome, optionally considering reverse complement and mismatches.
    
    Args:
        reads (list of str): Reads to match.
        genome (str): Reference genome.
        use_reverse_complement (bool): Whether to consider reverse complements.
        max_mismatches (int): Number of allowed mismatches.
    
    Returns:
        matched_count (int): Number of reads that matched at least once.
        total_reads (int): Total reads attempted.
    """
    matched_count = 0
    for read in reads:
        matches = []
        if max_mismatches == 0:
            if use_reverse_complement:
                matches = naive_search_with_reverse_complement(read, genome)
            else:
                matches = naive_exact_search(read, genome)
        else:
            # For approximate matching, currently no reverse complement check
            matches = naive_approx_search(read, genome, max_mismatches)
            if use_reverse_complement:
                rc = reverse_complement(read)
                matches += naive_approx_search(rc, genome, max_mismatches)
                matches = list(set(matches))
        if len(matches) > 0:
            matched_count += 1
    return matched_count, len(reads)

if __name__ == "__main__":
    genome_file = 'lambda_virus.fa'
    genome = read_genome(genome_file)

    # Generate artificial reads and test exact matching
    artificial_reads = generate_random_reads(genome, num_reads=100, read_len=100)
    matched, total = match_reads_to_genome(artificial_reads, genome)
    print(f"Artificial reads matched exactly: {matched} / {total}")

    # Load real reads from FASTQ for testing
    fastq_file = 'ERR266411_1.first1000.fastq'
    def read_fastq(filename):
        sequences = []
        with open(filename, 'r') as fh:
            while True:
                fh.readline()  # Skip header
                seq = fh.readline().rstrip()
                fh.readline()  # Skip plus line
                fh.readline()  # Skip quality line
                if len(seq) == 0:
                    break
                sequences.append(seq)
        return sequences

    real_reads = read_fastq(fastq_file)

    # Match real reads exactly
    matched, total = match_reads_to_genome(real_reads, genome)
    print(f"Real reads matched exactly: {matched} / {total}")

    # Match first 30 bases of reads exactly
    trimmed_reads = [r[:30] for r in real_reads]
    matched, total = match_reads_to_genome(trimmed_reads, genome)
    print(f"Real reads (first 30 bases) matched exactly: {matched} / {total}")

    # Match first 30 bases including reverse complements
    matched, total = match_reads_to_genome(trimmed_reads, genome, use_reverse_complement=True)
    print(f"Real reads (first 30 bases) matched exactly including reverse complement: {matched} / {total}")

    # Test approximate matching with mismatches
    pattern = 'TTCAAGCC'
    threshold = 3
    approx_matches = naive_approx_search(pattern, genome, threshold)
    print(f"Approximate matches for pattern '{pattern}' with up to {threshold} mismatches: {len(approx_matches)}")
    print(f"Positions: {approx_matches}")

    # Test naive search with reverse complement for a short pattern
    pattern2 = 'AGGT'
    matches_rc = naive_search_with_reverse_complement(pattern2, genome)
    print(f"Matches for pattern '{pattern2}' or reverse complement: {matches_rc}")
    print(f"Number of matches: {len(matches_rc)}")
