import pandas as pd
from math import log

# ============================================================
# 1. Known exon–intron boundary motif sequences (length = 9)
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

print("Motifs:")
for m in motifs:
    print(m)
print()

# ============================================================
# 2. Count Matrix
# ============================================================

count_matrix = {base: [0] * motif_length for base in alphabet}

for motif in motifs:
    for i, base in enumerate(motif):
        count_matrix[base][i] += 1

count_df = pd.DataFrame(count_matrix)
print("COUNT MATRIX")
print(count_df)
print()

# ============================================================
# 3. Relative Frequency (Weight) Matrix
# ============================================================

freq_matrix = {
    base: [count / num_motifs for count in counts]
    for base, counts in count_matrix.items()
}

freq_df = pd.DataFrame(freq_matrix)
print("RELATIVE FREQUENCY (WEIGHT) MATRIX")
print(freq_df)
print()

# ============================================================
# 4. Log-Likelihood Matrix
#    Background (Null Model): P(A)=P(C)=P(G)=P(T)=0.25
# ============================================================

background_prob = 0.25

log_likelihood_matrix = {
    base: [
        log(freq / background_prob) if freq > 0 else float("-inf")
        for freq in freqs
    ]
    for base, freqs in freq_matrix.items()
}

ll_df = pd.DataFrame(log_likelihood_matrix)
print("LOG-LIKELIHOOD MATRIX")
print(ll_df)
print()

# ============================================================
# 5. Analyze sequence S using sliding window
# ============================================================

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"
window_size = motif_length

def score_window(window):
    score = 0.0
    for i, base in enumerate(window):
        score += log_likelihood_matrix[base][i]
    return score

results = []

for i in range(len(S) - window_size + 1):
    window = S[i:i + window_size]
    score = score_window(window)
    results.append((i, window, score))

print("SLIDING WINDOW SCORES")
print("Position | Window     | Score")
print("----------------------------------")
for pos, win, sc in results:
    print(f"{pos:8d} | {win} | {sc:6.2f}")

# ============================================================
# 6. Best hit and biological interpretation
# ============================================================

best_hit = max(results, key=lambda x: x[2])

print("\nBEST SCORING WINDOW")
print(f"Position : {best_hit[0]}")
print(f"Sequence : {best_hit[1]}")
print(f"Score    : {best_hit[2]:.2f}")

if best_hit[2] > 0:
    print("\nINTERPRETATION:")
    print("✔ Strong signal detected.")
    print("✔ The sequence likely contains an exon–intron boundary.")
else:
    print("\nINTERPRETATION:")
    print("✘ No strong exon–intron boundary signal detected.")
