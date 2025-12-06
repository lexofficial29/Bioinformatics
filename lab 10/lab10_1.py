import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"

WINDOW = 30

def cg_content(seq):
    cg = seq.count("C") + seq.count("G")
    return (cg / len(seq)) * 100

def index_of_coincidence(seq):
    # IC = (sum(n_i * (n_i - 1))) / (N * (N - 1)) * 100
    counts = Counter(seq)
    N = len(seq)

    numerator = sum(n * (n - 1) for n in counts.values())
    denominator = N * (N - 1)
    return (numerator / denominator) * 100

CG_values = []
IC_values = []
centers = []

for i in range(len(S) - WINDOW + 1):
    window = S[i:i + WINDOW]

    CG_values.append(cg_content(window))
    IC_values.append(index_of_coincidence(window))

    # center of each pattern
    centers.append(i + WINDOW // 2)

print("CG for first window:", CG_values[0])   # should be ~29.27
print("IC for first window:", IC_values[0])   # should be ~27.53

plt.figure(figsize=(12, 5))
plt.plot(CG_values, label="C+G %")
plt.plot(IC_values, label="Index of Coincidence")
plt.title("DNA Pattern (Sliding Window = 30bp)")
plt.xlabel("Window position")
plt.ylabel("Value (%)")
plt.legend()
plt.grid(True)
plt.show()

# center of weight = weighted average of values based on position
# Use magnitude of pattern = CG + IC for weight
weights = np.array(CG_values) + np.array(IC_values)
positions = np.arange(len(weights))

center_of_weight = np.sum(weights * positions) / np.sum(weights)
print("\nCenter of Weight of Pattern =", center_of_weight)

plt.figure(figsize=(12, 5))
plt.scatter(centers, CG_values, s=20, label="CG% Centers")
plt.scatter(centers, IC_values, s=20, label="IC Centers")
plt.title("Centers of the Patterns")
plt.xlabel("Sequence Position")
plt.ylabel("Value (%)")
plt.legend()
plt.grid(True)
plt.show()