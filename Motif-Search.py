import random
import numpy as np

def score_kmer(profile, kmer):
    score = 1.0
    for i, nucleotide in enumerate(kmer):
        score *= profile['ACGT'.index(nucleotide)][i]
    return score

def probability_distribution(profile, kmers):
    prob_dict = {kmer: score_kmer(profile, kmer) for kmer in kmers}
    total = sum(prob_dict.values())
    return {kmer: prob / total for kmer, prob in prob_dict.items()}

def sample_kmer(prob_dist):
    kmers = list(prob_dist.keys())
    weights = list(prob_dist.values())
    return random.choices(kmers, weights=weights, k=1)[0]

def form_profile(motifs):
    k = len(motifs[0])
    t = len(motifs)
    profile = {i: [0] * k for i in range(4)}
    for i in range(k):
        column = [motif[i] for motif in motifs]
        for j, nucleotide in enumerate('ACGT'):
            profile[j][i] = (column.count(nucleotide) + 1) / (t + 4)
    return profile

def score_motifs(motifs):
    k = len(motifs[0])
    counts = np.zeros((4, k), dtype=int)
    for i in range(k):
        column = [motif[i] for motif in motifs]
        for j, nucleotide in enumerate('ACGT'):
            counts[j][i] = column.count(nucleotide)
    return int(sum(np.sum(counts, axis=0) - np.max(counts, axis=0)))

def form_motifs_from_profile(profile, dna, k):
    motifs = []
    for seq in dna:
        best_kmer, best_score = '', float('-inf')
        for j in range(len(seq) - k + 1):
            kmer = seq[j:j + k]
            score = score_kmer(profile, kmer)
            if score > best_score:
                best_score = score
                best_kmer = kmer
        motifs.append(best_kmer)
    return motifs

def gibbs_sampler(dna, k, t, N):
    motifs = [dna[i][random.randint(0, len(dna[i]) - k):][:k] for i in range(t)]
    best_motifs = list(motifs)
    best_score = float('inf')

    for _ in range(N):
        i = random.randrange(t)
        subset = motifs[:i] + motifs[i+1:]
        profile = form_profile(subset)
        kmers = [dna[i][j:j + k] for j in range(len(dna[i]) - k + 1)]
        prob_dist = probability_distribution(profile, kmers)
        motifs[i] = sample_kmer(prob_dist)
        current_score = score_motifs(motifs)
        if current_score < best_score:
            best_motifs, best_score = list(motifs), current_score

    return best_motifs, best_score

def gibbs_sampling(dna, k, t, N, runs):
    best_motifs, best_score = gibbs_sampler(dna, k, t, N)
    for _ in range(runs - 1):
        motifs, score = gibbs_sampler(dna, k, t, N)
        if score < best_score:
            best_motifs, best_score = motifs, score
    return best_motifs, best_score

def randomized_motif_search(dna, k, t, restarts):
    best_motifs, best_score = [], float('inf')
    for _ in range(restarts):
        motifs = [dna[i][random.randint(0, len(dna[i]) - k):][:k] for i in range(t)]
        while True:
            profile = form_profile(motifs)
            new_motifs = form_motifs_from_profile(profile, dna, k)
            current_score = score_motifs(new_motifs)
            if current_score < best_score:
                best_motifs, best_score = new_motifs, current_score
                motifs = new_motifs
            else:
                break
    return best_motifs, best_score

def read_dna_sequences(filepath):
    with open(filepath, 'r') as file:
        return file.read().strip().split()

def main():
    # Simple test DNA dataset
    dna = [
        "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
        "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
        "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
        "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
        "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
    ]
    k = 8
    t = len(dna)
    N = 100
    runs = 20

    print("=== Gibbs Sampling ===")
    gibbs_result, gibbs_score = gibbs_sampling(dna, k, t, N, runs)
    print("Best Motifs:", gibbs_result)
    print("Score:", gibbs_score)

    print("\n=== Randomized Motif Search ===")
    rand_result, rand_score = randomized_motif_search(dna, k, t, runs)
    print("Best Motifs:", rand_result)
    print("Score:", rand_score)

    # Additional unit-style tests
    print("\n=== Testing score_kmer ===")
    test_profile = {
        0: [0.2, 0.2, 0.3],
        1: [0.3, 0.3, 0.2],
        2: [0.2, 0.2, 0.3],
        3: [0.3, 0.3, 0.2]
    }
    test_kmer = "ACG"
    converted_profile = {i: test_profile[i] for i in range(4)}
    print("Score of", test_kmer, ":", score_kmer(converted_profile, test_kmer))

    print("\n=== Testing score_motifs ===")
    motifs = ["ATG", "ACG", "AAG"]
    print("Score:", score_motifs(motifs))

    print("\n=== Testing form_profile ===")
    prof = form_profile(motifs)
    print("Profile matrix:")
    for i, row in prof.items():
        print("ACGT"[i], ":", row)

if __name__ == "__main__":
    main()
