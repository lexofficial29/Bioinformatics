# Download 10 FASTA files containing different seq of influenza.
# Plot a chart for each column containing the 20 most freq repetitions of patterns.
from collections import Counter
import matplotlib.pyplot as plt
import os

# Function to find all 3–6 base pattern repeats in a sequence
def find_all_pattern_repeats(seq):
    pattern_counts = Counter()
    seq = seq.upper()
    for size in range(3, 7):
        for i in range(len(seq) - size):
            pattern = seq[i:i + size]
            repeats = 0
            pos = i
            while seq[pos:pos + size] == pattern:
                repeats += 1
                pos += size
            if repeats > 1:
                full_repeat = seq[i:i + repeats * size]
                pattern_counts[full_repeat] += 1
    return pattern_counts

# Directory with 10 text files
txt_dir = "fasta_files"
txt_files = [f for f in os.listdir(txt_dir) if f.endswith(".txt")]

# Create a figure with 10 subplots (adjust rows/cols as needed)
fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(20, 20))  # 5x2 grid
axes = axes.flatten()  # flatten 2D array for easy iteration

for ax, txt_file in zip(axes, txt_files):
    with open(os.path.join(txt_dir, txt_file), "r") as f:
        seq = f.read().replace("\n", "").strip()
    
    pattern_counts = find_all_pattern_repeats(seq)
    most_common = pattern_counts.most_common(20)
    
    if most_common:
        patterns, counts = zip(*most_common)
        ax.bar(range(len(patterns)), counts)
        ax.set_xticks(range(len(patterns)))
        ax.set_xticklabels(patterns, rotation=90, fontsize=8)
        ax.set_title(txt_file, fontsize=10)
        ax.set_ylabel("Frequency")
        ax.set_xlabel("Pattern")

plt.tight_layout()
plt.show()