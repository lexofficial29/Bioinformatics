# Show the alignment between the genomes in the console
# Make a visual representation of the alignment in which the regions that are similar are shown as long proportional rectangular bars.

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
from itertools import combinations
import time

# CHANGE THIS TO YOUR REAL EMAIL
Entrez.email = "your_email@example.com"

accessions = [
    "AF033819", "K03455", "K02013", "AY713414", "AJ249235",
    "AF005494", "AB098330", "AF077336", "AY371158", "U88824"
]

print("Downloading 10 HIV-1 genomes...\n")
genomes = []
names = []

for acc in accessions:
    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    genomes.append(record.seq)
    names.append(acc)
    print(f"{acc} → {len(record)} bp")
    handle.close()
    time.sleep(0.3)

print("\n" + "═" * 80)
print("   ALL PAIRWISE ALIGNMENTS WITH VISUAL BARS")
print("═" * 80 + "\n")

def smart_block_align(seq1, seq2, block_size=800, overlap=150):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    aligned1 = ""
    aligned2 = ""
    pos = 0
    while pos < max(len(seq1), len(seq2)):
        s1 = seq1[max(0, pos-overlap):pos+block_size]
        s2 = seq2[max(0, pos-overlap):pos+block_size]
        if len(s1) < 50 or len(s2) < 50:
            break
        alignments = aligner.align(s1, s2)
        a1, a2 = alignments[0][0], alignments[0][1]
        if pos > 0:
            a1, a2 = a1[overlap:], a2[overlap:]
        aligned1 += a1
        aligned2 += a2
        pos += block_size - overlap

    # Pad to same length
    max_len = max(len(aligned1), len(aligned2))
    aligned1 = aligned1.ljust(max_len, '-')
    aligned2 = aligned2.ljust(max_len, '-')
    return aligned1, aligned2

def make_visual_bar(aligned1, aligned2, width=80):
    bar1 = ""
    bar2 = ""
    matchline = ""
    for a, b in zip(aligned1, aligned2):
        if a == b and a != '-':
            bar1 += "█"
            bar2 += "█"
            matchline += "|"
        elif a == '-' or b == '-':
            bar1 += "░"
            bar2 += "░"
            matchline += " "
        else:
            bar1 += "░"
            bar2 += "░"
            matchline += "·"

    step = max(1, len(bar1) // width)
    return bar1[::step][:width], bar2[::step][:width], matchline[::step][:width]

# MAIN LOOP — ALL 45 PAIRS
for (i, j) in combinations(range(10), 2):
    acc1, acc2 = names[i], names[j]
    seq1, seq2 = genomes[i], genomes[j]

    print(f"{acc1}  vs  {acc2}")
    
    aln1, aln2 = smart_block_align(seq1, seq2)

    # Calculate identity
    matches = sum(a == b and a != '-' for a, b in zip(aln1, aln2))
    total = sum(a != '-' for a in aln1)
    identity = 100 * matches / total if total > 0 else 0

    # CORRECTED: Build match line properly
    preview_len = 120
    match_line = ''.join('|' if aln1[k] == aln2[k] and aln1[k] != '-' else '.' 
                         for k in range(min(preview_len, len(aln1))))

    print("Alignment preview:")
    print(f"1: {aln1[:preview_len]}")
    print(f"   {match_line}")
    print(f"2: {aln2[:preview_len]}")
    print()

    # Visual bars
    bar1, bar2, vis_match = make_visual_bar(aln1, aln2, width=80)
    print("Visual similarity map:")
    print(f"1: {bar1}")
    print(f"   {vis_match}")
    print(f"2: {bar2}")
    print(f"   Identity: {identity:.2f}%\n")
    print("═" * 80 + "\n")

print("All 45 pairwise visual alignments completed!")