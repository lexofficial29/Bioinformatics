# Download 10 influenza genomes and use the restriction enzimes from the set to simulate the electrophoresis gel for these viral genomes.
# Eliminate all communalities from the 10 electrophoresis gel genomes and make 1 gemneral electrophoresis gel that shows only the differences between these genomes. 

import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from collections import Counter

folder = "fasta_files"

enzymes = [
    ("EcoRI",   "GAATTC", 1),
    ("BamHI",   "GGATCC", 1),
    ("HindIII", "AAGCTT", 1),
    ("TaqI",    "TCGA",   1),
    ("HaeIII",  "GGCC",   2),
]

strain_names = ["A/California", "A/Perth", "A/Vietnam", "A/duck", "A/Hong Kong",
                "A/Texas", "A/Yamagata", "A/Victoria", "A/New York", "A/chicken"]

def read_dna(i):
    with open(os.path.join(folder, f"{i}.txt")) as f:
        return "".join(line.strip() for line in f if not line.startswith(">")).upper()

def digest(seq, site, cut_after):
    cuts = [i + cut_after for i in range(len(seq)-len(site)+1) if seq[i:i+len(site)] == site]
    if not cuts:
        return [len(seq)]
    cuts = sorted(cuts)
    frags = [cuts[0]] + [b-a for a,b in zip(cuts,cuts[1:])] + [len(seq)-cuts[-1]]
    return sorted(frags, reverse=True)

def mig(size):
    return 9.8 - 7.8 * (np.log10(max(size,50)) - 1.7) / (np.log10(20000) - 1.7)

ladder = [10000,8000,6000,5000,4000,3000,2000,1500,1000,750,500,250,100]

all_combined = []

print("="*80)
print("RESTRICTION DIGEST RESULTS – 10 INFLUENZA GENOMES")
print("="*80)

for i in range(1, 11):
    seq = read_dna(i)
    print(f"\nGenome {i:2d} → {strain_names[i-1]:18} → {len(seq):,} bp")
    
    per_enzyme = []
    combined_set = set()

    for name, site, off in enzymes:
        frags = digest(seq, site, off)
        cuts = len(frags) - 1 if len(seq) in frags else len(frags)
        print(f"  {name:8} → {cuts:2} cut(s) → fragments: {frags}")
        per_enzyme.append(frags)
        combined_set.update(frags)

    combined_list = sorted(combined_set, reverse=True)
    all_combined.append(combined_set)

    # Plot individual gel
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_facecolor('black')

    # Ladder
    for s in ladder:
        y = mig(s)
        ax.add_patch(Rectangle((0.4, y-0.06), 0.25, 0.12, color='orange', alpha=0.9))
    ax.text(0.525, 9.7, "1 kb\nLadder", color='orange', ha='center', weight='bold', fontsize=11)

    # Lanes
    lane_labels = [e[0] for e in enzymes] + ["Combined"]
    lane_data = per_enzyme + [combined_list]

    for lane_idx, frags in enumerate(lane_data):
        x = 1.8 + lane_idx * 1.15
        for size in frags:
            if size < 50: continue
            y = mig(size)
            thick = max(0.07, np.sqrt(size)/38)
            color = 'cyan' if lane_idx < 5 else 'lime'
            ax.add_patch(Rectangle((x-0.4, y-thick/2), 0.8, thick,
                                  facecolor=color, edgecolor='white', lw=0.9))

    ax.set_xlim(0, 10.5)
    ax.set_ylim(0, 10)
    ax.set_xticks([1.8 + i*1.15 for i in range(6)])
    ax.set_xticklabels(lane_labels, color='white', fontsize=12)
    ax.set_yticks([])
    ax.set_title(f"Strain {i} — {strain_names[i-1]}\n5-enzyme restriction digest", 
                 color='white', fontsize=14, pad=20)

    plt.tight_layout()
    plt.savefig(f"strain_{i}.png", dpi=300, facecolor='black', bbox_inches='tight')
    plt.close()
    print(f"  → strain_{i}.png saved\n")

# Final differences summary + image
print("="*80)
print("DIFFERENCES BETWEEN THE 10 STRAINS")
print("="*80)

all_sizes = [s for sett in all_combined for s in sett]
counter = Counter(all_sizes)
diff = sorted([s for s,c in counter.items() if c < 10], reverse=True)

if not diff:
    print("All 10 strains have IDENTICAL restriction patterns!")
else:
    for rank, size in enumerate(diff, 1):
        print(f"{rank:2d}. {size:5} bp  ← present in {counter[size]:2}/10 strains")

# Save differences-only gel
fig, ax = plt.subplots(figsize=(7, 9), facecolor='black')
ax.set_facecolor('black')

for s in ladder:
    y = mig(s)
    ax.add_patch(Rectangle((0.4, y-0.06), 0.25, 0.12, color='orange'))
ax.text(0.525, 9.7, "1 kb\nLadder", color='orange', ha='center', weight='bold')

x=12

for size in diff:
    if size < 50: continue
    y = mig(size)
    thick = max(0.12, np.sqrt(size)/28)
    ax.add_patch(Rectangle((2-0.5, y-thick/2), 1.0, thick,
                          facecolor='lime', edgecolor='white', lw=2))

ax.set_xlim(0, 4)
ax.set_ylim(0, 10)
ax.set_xticks([2])
ax.set_xticklabels(["Variable\nbands"], color='lime', fontsize=14)
ax.set_yticks([])
ax.set_title("Only fragments that differ among the 10 influenza strains", color='white', fontsize=15)

plt.tight_layout()
plt.savefig("differences_only.png", dpi=300, facecolor='black', bbox_inches='tight')
plt.close()