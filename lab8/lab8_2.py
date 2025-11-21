import random

genome = "AGGAGAATAGCCACTTCTGATACTGGAGAGCATGGGTCGTCGATTAGCTTTGTCCATAGAGCTGTCCCACTCTGTCGCATAAAAGCTGAAGCGTTTTAACCGGATTAACCGGCCAGTCCGCATGCGTACCTAGGTATTCGAAGGCCCGGGAAATTTTTCGGGATTTACCGCGTACCCCCGGGAAATTTTGAGGGATTTACCGCATGTATAGTTAACCGGATTAACCCCATAATCTGTCGCTACGCCGATC"

transposable = {
    "TE1": "ATGCGTACCTAGGTA",
    "TE2": "GGGATTTACCG",
    "TE3": "TTAACCGGATTAACC",
    "TE4": "CCCGGGAAATTT" 
}

def find_known_TEs(genome, transposable):
    results = []
    for name, seq in transposable.items():
        start = genome.find(seq)
        while start != -1:
            results.append((name, start, start+len(seq)))
            start = genome.find(seq, start + 1)
    return results

hits = find_known_TEs(genome, transposable)

print("Detected Transposable Elements in the genome:")
for name, start, end in hits:
    print(f"{name}: start={start}, end={end}, length={end-start}")
