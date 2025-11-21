import random

length = random.randint(200, 400)
bases = ["A", "C", "G", "T"]
genome = [random.choice(bases) for _ in range(length)]

transposable = {
    "TE1": "ATGCGTACCTAGGTA",
    "TE2": "GGGATTTACCG",
    "TE3": "TTAACCGGATTAACC",
    "TE4": "CCCGGGAAATTT" 
}

inserted_positions = {}

for name, seq in transposable.items():
    for _ in range(2):
        pos = random.randint(0, len(genome) - len(seq))
        if name not in inserted_positions:
            inserted_positions[name] = []
        inserted_positions[name].append((pos, pos + len(seq)))

        genome[pos : pos + len(seq)] = list(seq)

genome = "".join(genome)

print("Genome length:", len(genome))
print("\nArtificial Genome Sequence\n")
print(genome)

print("\nInserted Transposable Element Positions")
for name, positions in inserted_positions.items():
    for start, end in positions:
        print(f"{name}: start={start}, end={end}, length={end-start}")
