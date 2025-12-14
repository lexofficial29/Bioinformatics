# Downland 10 Hiv genomes and use them as an input for your application. Adapt your current algorithm of peers wise alignment of 2 genomes. 
# Note: take into consideration the fact that your current implementation cannot align 2 games directly, therefore, you must find a strategy im which the 2 genoms can be aligned step by step on smoller regions. Make the alligment between each genome

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
from itertools import combinations
import os
import time

# 10 diverse HIV-1 complete genomes (from different subtypes and countries)
accessions = [
    "AF033819",   # Subtype B (USA)
    "K03455",     # HXB2 reference (Subtype B)
    "K02013",     # LAV-1/BRU (Subtype B)
    "AY713414",   # Subtype C (South Africa)
    "AJ249235",   # Subtype A (Kenya)
    "AF005494",   # Subtype D (Uganda)
    "AB098330",   # Subtype CRF01_AE (Thailand)
    "AF077336",   # Subtype G (Russia)
    "AY371158",   # Subtype F1 (Brazil)
    "U88824",     # Subtype O (Cameroon outlier group)
]

os.makedirs("hiv_pairwise_alignments", exist_ok=True)

print("Downloading 10 HIV-1 complete genomes from NCBI...")
genomes = []
names = []

for acc in accessions:
    print(f"  Downloading {acc}...", end="")
    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    record.id = acc
    record.description = ""
    genomes.append(record)
    names.append(acc)
    print(f" {len(record)} bp")
    handle.close()
    time.sleep(0.5)  # Be nice to NCBI

print(f"\nDownloaded {len(genomes)} genomes.\n")

def block_pairwise_align(seq1, seq2, block_size=600, overlap=100):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -7
    aligner.extend_gap_score = -2

    len1, len2 = len(seq1), len(seq2)
    aligned1_parts = []
    aligned2_parts = []

    pos = 0
    while pos < max(len1, len2):
        start1 = max(0, pos - overlap)
        end1 = min(len1, pos + block_size)
        start2 = max(0, pos - overlap)
        end2 = min(len2, pos + block_size)

        block1 = seq1[start1:end1]
        block2 = seq2[start2:end2]

        if len(block1) < 10 or len(block2) < 10:
            break

        alignments = aligner.align(block1, block2)
        best = alignments[0]
        al1, al2 = str(best).split('\n')[:2]

        if pos > 0:
            al1 = al1[overlap:]
            al2 = al2[overlap:]

        aligned1_parts.append(al1)
        aligned2_parts.append(al2)

        pos += (block_size - overlap)

    final_al1 = ''.join(aligned1_parts)
    final_al2 = ''.join(aligned2_parts)

    diff = abs(len(final_al1) - len(final_al2))
    if len(final_al1) < len(final_al2):
        final_al1 += '-' * diff
    elif len(final_al2) < len(final_al1):
        final_al2 += '-' * diff

    return final_al1, final_al2

print("Starting pairwise alignment of all 45 genome pairs...\n")
print(f"{'Pair':<25} {'Identity (%)':<12} {'Aligned Length':<15}")
print("-" * 60)

results = []

for i, j in combinations(range(len(genomes)), 2):
    acc1, acc2 = names[i], names[j]
    seq1, seq2 = genomes[i].seq, genomes[j].seq

    print(f"Aligning {acc1} vs {acc2} ... ", end="")

    aligned1, aligned2 = block_pairwise_align(seq1, seq2, block_size=600, overlap=100)

    # Calculate identity
    matches = sum(a == b for a, b in zip(aligned1, aligned2) if a != '-' and b != '-')
    total_aligned = len(aligned1.replace('-', ''))
    identity = (matches / total_aligned) * 100 if total_aligned > 0 else 0

    print(f"{identity:6.2f}%")

    # Save alignment to file
    filename = f"hiv_pairwise_alignments/{acc1}_vs_{acc2}.aln.txt"
    with open(filename, "w") as f:
        f.write(f">{acc1}\n{aligned1}\n\n>{acc2}\n{aligned2}\n")
        f.write(f"\nIdentity: {identity:.2f}% over {total_aligned} bases\n")

    results.append((acc1, acc2, identity, total_aligned))

print("\nAll 45 pairwise alignments completed!")
print("Files saved in: hiv_pairwise_alignments/\n")

print("Summary (sorted by sequence identity):")
print(f"{'Rank':<4} {'Genome 1':<12} {'Genome 2':<12} {'Identity (%)'}")
print("-" * 50)
for rank, (g1, g2, ident, length) in enumerate(sorted(results, key=lambda x: x[2], reverse=True), 1):
    print(f"{rank:<4} {g1:<12} {g2:<12} {ident:8.2f}%")