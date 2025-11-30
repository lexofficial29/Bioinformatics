# 1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).
# 2. Use 5 restriction enzymes (enzyme name, recognized sequence, cleavage site):
# Input of the implementation:
# 1. The recognized sequence for each restriction enzyme.
# 2. A DNA sequence to be digested.

# Output of the implementation:
# 1. Number of cleavages for a each restriction enzyme.
# 2. Cleavage positions and length of fragments.
# 3. A simulation of the electrophoresis gel based on the number of restriction enzymes used.

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import random

dna = "ATGGCGTCTCAAGGCACCAAACGATCTTATGAACAGATGGAAACTGGTGGAGAACGCCAGAATGCCACTGAAATCAGAGCATCTGTTGGAAGAATGGTTGGTGGAATTGGAAGGTTTTATATACAGATGTGCACTGAACTCAAACTCAGCGATTATGAGGGGAGACTGATCCAGAACAGCATAACAATAGAGAGAATGATTCTCTCTGCATTTGATGAAAGGAGGAACAAGTACCTGGAAGAACATCCCAGTGCGGGAAAGGACCCAAAGAAAACTGGAGGTCCAATCTACAGAAGAAGAGACGGAAAGTGGATGAGGGAGCTGATTCTGTATGACAAAGAAGAGATCAGAAGGATCTGGCGTCAAGCGAATAATGGAGAAGATGCAACTGCTGGTCTCACCCATCTGATGATCTGGCACTCCAACCTGAATGATGCCACATATCAGAGAACAAGGGCTCTCGTGCGCACTGGAATGGATCCCAGAATGTGCTCTCTGATGCAAGGATCAACTCTCCCAAGAAGGTCTGGAGCTGCTGGTGCAGCAGTAAAAGGGGTCGGAACAATGGTAATGGAATTGATTAGAATGATAAAGCGAGGGATTAATGATCGGAATTTCTGGAGAGGCGAAAATGGAAGAAGGACAAGGATTGCCTATGAGAGAATGTGCAACATCCTCAAAGGGAAATTTCAAACAGCAGCACAAAGAGCAATGATGGATCAAGTGCGAGAAAGCAGGAATCCTGGGAATGCTGAAATTGAAGATCTCATTTTTCTGGCACGGTCTGCACTCATCCTGAGAGGGTCAGTGGCCCACAAGTCTTGTCTGCCTGCTTGTGTTTACGGACTTGCTGTGGCCAGTGGATACGACTTTGAGAGAGAAGGATACTCTCTGGTCGGAATAGACCCTTTCCGTCTGCTTCAAAACAGCCAGGTCTTCAGTCTCATTAGACCAAATGAAAACCCAGCACATAAAAGCCAATTGGTATGGATGGCATGCCATTCAGCAGCGTTTGAGGACCTGAGGGTATCAAGTTTCATCAGAGGGACAAAAGTGGTCCCAAGAGGACAACTATCCACCAGAGGAGTTCAAATTGCATCAAATGAAAACATGGAAACAATGGACTCCAGCACTCTTGAATTGAGAAGCAGATACTGGGCTATAAGAACCAGGAGTGGAGGAAACACCAACCAACAGAGAGCTTCTGCAGGACAAATCAGCGTACAACCCACCTTCTCAGTACAGAGAAATCTTCCCTTTGAAAGAGCGACCATCATGGCGGCATTTACAGGGAACACTGAAGGCAGGACCTCTGACATGAGGACTGAAATCATAAGAATGATGGAAAGTGCCAAACCAGAAGATGTGTCCTTCCAGGGGCGGGGAGTCTTCGAGCTCTCGGACGAAAAGGCAACGAACCCGATCGTGCCTTCCTTTGACATGAGCAACGAAGGATCTTATTTCTTCGGAGACAGTGCAGAGGAGTATGACAATTAA"

enzymes = [
    ("EcoRI",   "GAATTC", 1),   # cuts after 1st base → G^AATTC
    ("BamHI",   "GGATCC", 1),   # G^GATCC
    ("HindIII", "AAGCTT", 1),   # A^AGCTT
    ("TaqI",    "TCGA",   1),   # T^CGA
    ("HaeIII",  "GGCC",   2),   # GG^CC
]

def find_cut_sites(seq, site, cut_after):
    positions = []
    n = len(site)
    for i in range(len(seq) - n + 1):
        if seq[i:i+n] == site:
            positions.append(i + cut_after)  # 0-based index + offset → cut position
    return positions

def get_fragments(length, cuts):
    if not cuts:
        return [length]
    cuts = sorted(cuts)
    fragments = [cuts[0]] + [cuts[i] - cuts[i-1] for i in range(1, len(cuts))] + [length - cuts[-1]]
    return sorted(fragments, reverse=True)

results = []
print("RESTRICTION DIGEST RESULTS")
print("="*60)
for name, site, offset in enzymes:
    cuts = find_cut_sites(dna, site, offset)
    frags = get_fragments(len(dna), cuts)
    results.append((name, frags))
    print(f"{name:8} → {len(cuts)} cut(s) → fragments: {frags}")

# Add uncut lane
results.insert(0, ("Uncut", [len(dna)]))

plt.style.use('dark_background')
fig, ax = plt.subplots(figsize=(11, 8))

# Log-like migration function
def migration_distance(size_bp):
    return 9.5 - 7.5 * (np.log10(size_bp) - 1) / (np.log10(10000) - 1)

# 1 kb ladder bands
ladder_sizes = [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1200, 1000,
                800, 600, 500, 400, 300, 200, 100, 50]

lane_positions = np.arange(len(results))
lane_width = 0.7

# Draw each lane
for lane_idx, (label, fragments) in enumerate(results):
    x = lane_positions[lane_idx]
    
    # Lane background
    ax.add_patch(Rectangle((x - lane_width/2, 0), lane_width, 10,
                          color='#0a0a1a', alpha=0.8, zorder=1))
    
    # Draw DNA bands
    for size in fragments:
        if size < 40:  # too small to see
            continue
        y = migration_distance(size)
        thickness = max(0.08, np.sqrt(size) / 40)
        brightness = min(1.0, size / 2000)
        
        ax.add_patch(Rectangle((x - lane_width*0.8/2, y - thickness/2),
                               lane_width*0.8, thickness,
                               facecolor='cyan', edgecolor='white',
                               linewidth=0.8, alpha=brightness, zorder=3))

# Draw 1 kb DNA ladder (left side)
ax.text(-1.3, 9.7, "1 kb\nLadder", color="orange", fontsize=12, ha='center', weight='bold')
for size in ladder_sizes:
    y = migration_distance(size)
    ax.add_patch(Rectangle((-1.6, y - 0.08), 0.8, 0.16,
                           facecolor='orange', edgecolor='yellow',
                           linewidth=1, alpha=0.9, zorder=4))
    if size >= 1000:
        ax.text(-2.2, y, f"{size//1000} kb", color="orange", va='center', fontsize=9)
    else:
        ax.text(-2.2, y, f"{size}", color="orange", va='center', fontsize=9)

# Styling
ax.set_xlim(-3, len(results))
ax.set_ylim(0, 10)
ax.set_xticks(lane_positions)
ax.set_xticklabels([name for name, _ in results], rotation=40, ha='right', color='white', fontsize=11)
ax.set_yticks([])
ax.set_title("Restriction Enzyme Digest – Agarose Gel Simulation", 
             color='white', fontsize=18, pad=30)
ax.spines[['top','right','left','bottom']].set_visible(False)

plt.tight_layout()
plt.show()