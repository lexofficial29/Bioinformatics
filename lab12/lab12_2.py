# Download 10 influenza genomes (i have them in a folder named fasta_files named from 1 to 10.txt) 
# Adapt your application from the previous assignment in order to scan each genome for possible motives.
# For each genome make a chart that shows the signal with the most likely locations of real functional motifs

import os
import matplotlib.pyplot as plt
import pandas as pd
from math import log

# ============================================================
# 1. Motifs reconstructed COLUMN-WISE from previous assignment
# ============================================================

motifs = [
    "GTCATTACTA",
    "ACACAATAGA",
    "GCGAGGGGTG",
    "GGGGGGGGGG",
    "TTTTTTTTTT",
    "AATCCAAAGA",
    "AAGAACATAA",
    "AGGGTTCAGG",
    "CTATTGTCTT"
]

alphabet = ["A", "C", "G", "T"]
motif_length = len(motifs[0])
num_motifs = len(motifs)

# ============================================================
# 2. Build Count Matrix
# ============================================================

count_matrix = {base: [0] * motif_length for base in alphabet}

for motif in motifs:
    for i, base in enumerate(motif):
        count_matrix[base][i] += 1

# ============================================================
# 3. Relative Frequency Matrix
# ============================================================

freq_matrix = {
    base: [count / num_motifs for count in counts]
    for base, counts in count_matrix.items()
}

# ============================================================
# 4. Log-Likelihood Matrix (background = 0.25)
# ============================================================

background = 0.25

log_likelihood = {
    base: [
        log(freq / background) if freq > 0 else float("-inf")
        for freq in freqs
    ]
    for base, freqs in freq_matrix.items()
}

# ============================================================
# 5. Sliding window scoring function
# ============================================================

def score_window(window):
    score = 0.0
    for i, base in enumerate(window):
        if base in log_likelihood:
            score += log_likelihood[base][i]
        else:
            return float("-inf")
    return score

# ============================================================
# 6. Read genome from FASTA or raw text
# ============================================================

def read_genome(filepath):
    seq = ""
    with open(filepath, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip().upper()
    return seq

# ============================================================
# 7. Scan genome and return scores
# ============================================================

def scan_genome(sequence):
    scores = []
    positions = []

    for i in range(len(sequence) - motif_length + 1):
        window = sequence[i:i + motif_length]
        score = score_window(window)
        scores.append(score)
        positions.append(i)

    return positions, scores

# ============================================================
# 8. Process all 10 influenza genomes
# ============================================================

input_folder = "fasta_files"
output_folder = "plots"

os.makedirs(output_folder, exist_ok=True)

for i in range(1, 11):
    file_path = os.path.join(input_folder, f"{i}.txt")
    genome = read_genome(file_path)

    positions, scores = scan_genome(genome)

    max_score = max(scores)
    threshold = max_score * 0.8  # heuristic threshold

    # ========================================================
    # 9. Plot signal
    # ========================================================

    plt.figure(figsize=(12, 4))
    plt.plot(positions, scores)
    plt.axhline(threshold, linestyle="--")

    plt.title(f"Influenza Genome {i} – Motif Signal")
    plt.xlabel("Genome Position")
    plt.ylabel("Log-Likelihood Score")

    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"genome_{i}_signal.png"))
    plt.close()

    print(f"Genome {i}: max score = {max_score:.2f}")

print("\nAnalysis complete. Signal plots saved in 'plots/' folder.")