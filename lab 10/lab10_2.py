# Download 10 influenza genomes and 10 covid 19 gemone. use their fasta files to plot their objective digital straints.
# on a second type of chart plot thecenter of weight for each ODS. label each of the 20 centers of weight(plotted as a small circle) with the name and variant of the viruses.

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from Bio import Entrez, SeqIO
import time
import os

Entrez.email = "alyoblu@gmail.com"
WINDOW = 300
OUTPUT_DIR = "viral_genomes"
os.makedirs(OUTPUT_DIR, exist_ok=True)

genomes = [
    # === Influenza A (various subtypes and years) ===
    {"accession": "CY077311", "name": "A/Puerto Rico/8/1934 (H1N1)"},
    {"accession": "CY021853", "name": "A/California/04/2009 (H1N1pdm09)"},
    {"accession": "MK624330", "name": "A/Hong Kong/4801/2014 (H3N2)"},
    {"accession": "KX867461", "name": "A/mallard/Ohio/16OS1072/2016 (H5N8)"},
    {"accession": "MF194110", "name": "A/Vietnam/1203/2004 (H5N1)"},
    {"accession": "CY009348", "name": "A/Brevig Mission/1/1918 (H1N1)"},
    {"accession": "LC333185", "name": "A/duck/Mongolia/129/2015 (H5N1)"},
    {"accession": "MT419816", "name": "A/Shanghai/2/2013 (H7N9)"},
    {"accession": "KJ942682", "name": "A/chicken/Jiangsu/1022/2013 (H10N8)"},
    {"accession": "CY166976", "name": "A/New York/392/2004 (H3N2)"},
    
    # === SARS-CoV-2 (COVID-19) - key variants ===
    {"accession": "NC_045512", "name": "Wuhan-Hu-1 (Original)"},
    {"accession": "OM212469", "name": "Alpha (B.1.1.7)"},
    {"accession": "OU296679", "name": "Beta (B.1.351)"},
    {"accession": "OV108635", "name": "Gamma (P.1)"},
    {"accession": "OK104651", "name": "Delta (B.1.617.2)"},
    {"accession": "OM617922", "name": "Omicron BA.1"},
    {"accession": "OP927946", "name": "Omicron BA.2"},
    {"accession": "OR164349", "name": "Omicron BA.5"},
    {"accession": "OX382250", "name": "Omicron XBB.1.5"},
    {"accession": "PP298089", "name": "Omicron JN.1 (2024)"},
]

def download_genome(accession, filename):
    """Download genome from NCBI if not already downloaded"""
    if os.path.exists(filename):
        return filename
    print(f"Downloading {accession}...")
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        with open(filename, "w") as f:
            f.write(handle.read())
        time.sleep(1)  # Be nice to NCBI
        return filename
    except Exception as e:
        print(f"Failed to download {accession}: {e}")
        return None

def cg_content(seq):
    cg = seq.count("C") + seq.count("G")
    return (cg / len(seq)) * 100 if len(seq) > 0 else 0

def index_of_coincidence(seq):
    if len(seq) < 2:
        return 0
    counts = Counter(seq.upper())
    N = len(seq)
    numerator = sum(n * (n - 1) for n in counts.values())
    return (numerator / (N * (N - 1))) * 100

def compute_ods(sequence, window=300):
    seq = str(sequence).upper().replace("-", "").replace("N", "A")  # Clean
    n = len(seq)
    if n < window:
        return [], []
    
    cg_vals = []
    ic_vals = []
    centers = []
    
    for i in range(n - window + 1):
        win = seq[i:i+window]
        cg_vals.append(cg_content(win))
        ic_vals.append(index_of_coincidence(win))
        centers.append(i + window // 2)
    
    ods = np.array(cg_vals) + np.array(ic_vals)
    return centers, ods

results = []

plt.figure(figsize=(14, 8))

for idx, genome in enumerate(genomes):
    acc = genome["accession"]
    name = genome["name"]
    filename = os.path.join(OUTPUT_DIR, f"{acc}.fasta")
    
    path = download_genome(acc, filename)
    
    if not path:
        print(f"Skipping {acc}")
        continue
    
    record = next(SeqIO.parse(path, "fasta"))
    seq = str(record.seq)
    
    centers, ods = compute_ods(seq, WINDOW)
    
    # Store for center of weight
    if len(ods) > 0:
        weights = ods
        positions = np.array(centers)
        center_of_weight = np.sum(weights * positions) / np.sum(weights)
        
        results.append({
            "name": name,
            "accession": acc,
            "center": center_of_weight,
            "length": len(seq),
            "ods": ods,
            "positions": centers
        })
        
        # Plot ODS curve
        color = "blue" if "A/" in name else "red"
        alpha = 0.6 if "A/" in name else 0.8
        plt.plot(centers, ods, label=name if idx < 5 else None, color=color, alpha=alpha, linewidth=1)

plt.title("Objective Digital Strain (ODS = CG% + IC) - Sliding Window 300 bp", fontsize=14)
plt.xlabel("Genomic Position (nt)")
plt.ylabel("ODS Value (%)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 8))
flu_centers = [r["center"] for r in results if "A/" in r["name"]]
flu_names = [r["name"] for r in results if "A/" in r["name"]]

cov_centers = [r["center"] for r in results if "SARS" in r["name"] or "Wuhan" in r["name"] or "Omicron" in r["name"]]
cov_names = [r["name"] for r in results if "SARS" in r["name"] or "Wuhan" in r["name"] or "Omicron" in r["name"]]

plt.scatter(flu_centers, [1]*len(flu_centers), s=100, c="blue", label="Influenza A", alpha=0.8)
plt.scatter(cov_centers, [2]*len(cov_centers), s=120, c="red", label="SARS-CoV-2", alpha=0.9)

# Annotate each point
for r in results:
    x = r["center"]
    y = 1 if "A/" in r["name"] else 2
    plt.annotate(r["name"].split("/")[-1] if "/" in r["name"] else r["name"][:20],
                 (x, y),
                 xytext=(5, 10), textcoords='offset points',
                 fontsize=9, ha='left', rotation=0)

plt.yticks([1, 2], ["Influenza A", "SARS-CoV-2"])
plt.xlabel("Center of Weight (position in nt)")
plt.title("Center of Weight of Objective Digital Strain (ODS)\nfor 10 Influenza A and 10 SARS-CoV-2 Genomes")
plt.grid(True, axis='x', alpha=0.3)
plt.tight_layout()
plt.show()

# Print summary
print("\n=== Center of Weight Summary ===")
for r in results:
    print(f"{r['name']:45} → {r['center']:.0f} nt (length: {r['length']} nt)")