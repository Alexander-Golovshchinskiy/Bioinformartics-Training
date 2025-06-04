import matplotlib.pyplot as plt
from collections import Counter

def read_fastq(filename):
    """
    Read sequences and quality strings from a FASTQ file.
    
    Args:
        filename (str): Path to the FASTQ file.
    
    Returns:
        sequences (list of str): List of DNA sequences.
        qualities (list of str): List of quality strings (Phred+33 encoded).
    """
    sequences = []
    qualities = []
    with open(filename, 'r') as fh:
        while True:
            fh.readline()  # Header line, ignored
            seq = fh.readline().rstrip()
            fh.readline()  # Plus line, ignored
            qual = fh.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def phred33_to_q(char):
    """
    Convert a single Phred+33 character to quality score (integer).
    """
    return ord(char) - 33

def qualities_to_scores(quality_list):
    """
    Convert a list of Phred+33 quality strings to numeric quality scores.
    
    Args:
        quality_list (list of str): List of quality strings.
    
    Returns:
        scores (list of list of int): Numeric quality scores.
    """
    all_scores = []
    for qual_str in quality_list:
        scores = [phred33_to_q(ch) for ch in qual_str]
        all_scores.append(scores)
    return all_scores

def create_quality_histogram(quality_list, max_score=50):
    """
    Create a histogram of quality scores from a list of Phred+33 quality strings.
    
    Args:
        quality_list (list of str): List of quality strings.
        max_score (int): Max expected quality score (default 50).
    
    Returns:
        hist (list of int): Counts of each quality score index.
    """
    hist = [0] * max_score
    for qual_str in quality_list:
        for ch in qual_str:
            q = phred33_to_q(ch)
            if q < max_score:
                hist[q] += 1
    return hist

def calculate_gc_content_by_position(sequences):
    """
    Calculate GC content for each base position across a list of sequences.
    
    Args:
        sequences (list of str): List of equal-length DNA sequences.
    
    Returns:
        gc_content (list of float): GC proportion at each position.
    """
    seq_len = len(sequences[0])
    gc_counts = [0] * seq_len
    totals = [0] * seq_len
    
    for seq in sequences:
        for i, base in enumerate(seq):
            if base in ('G', 'C'):
                gc_counts[i] += 1
            totals[i] += 1
    
    gc_content = [gc_counts[i] / totals[i] if totals[i] > 0 else 0 for i in range(seq_len)]
    return gc_content

def plot_gc_content(gc_content, ax=None):
    """
    Plot GC content by position.
    
    Args:
        gc_content (list of float): GC content per position.
        ax (matplotlib axis, optional): Axis to plot on.
    """
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(range(len(gc_content)), gc_content)
    ax.set_xlabel('Position in read')
    ax.set_ylabel('GC content')
    ax.set_title('GC content by position')
    plt.tight_layout()

def plot_quality_histogram(hist):
    """
    Plot histogram of quality scores.
    
    Args:
        hist (list of int): Histogram counts of quality scores.
    """
    plt.bar(range(len(hist)), hist)
    plt.xlabel('Quality score (Phred33)')
    plt.ylabel('Count')
    plt.title('Quality score distribution')
    plt.show()

def count_bases(sequences):
    """
    Count base frequencies across all sequences.
    
    Args:
        sequences (list of str): List of DNA sequences.
    
    Returns:
        Counter object with base counts.
    """
    count = Counter()
    for seq in sequences:
        count.update(seq)
    return count

if __name__ == "__main__":
    # Example usage
    fastq_file = 'ERR037900_1.first1000.fastq'
    sequences, qualities = read_fastq(fastq_file)

    print(f"Number of sequences: {len(sequences)}")
    print(f"Length of first sequence: {len(sequences[0])}")
    print(f"First 5 quality strings: {qualities[:5]}")

    # Convert qualities to scores (not printed fully for clarity)
    quality_scores = qualities_to_scores(qualities)

    # Create and plot quality histogram
    quality_hist = create_quality_histogram(qualities)
    plot_quality_histogram(quality_hist)

    # Calculate and plot GC content
    gc_content = calculate_gc_content_by_position(sequences)
    plot_gc_content(gc_content)

    # Count base frequencies and print
    base_counts = count_bases(sequences)
    print("Base counts across all sequences:", base_counts)
