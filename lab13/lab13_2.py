# 1.Use a random DNA seq of about 50 letters. Use the seq to compute the transition probabilities between letters.
# Your output should be the transition matrix stored as a JSON file.

import random
import json
from collections import defaultdict

def generate_dna_sequence(length=50):
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))


def compute_transition_matrix(sequence):
    transitions = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)

    for i in range(len(sequence) - 1):
        current_base = sequence[i]
        next_base = sequence[i + 1]
        transitions[current_base][next_base] += 1
        totals[current_base] += 1

    transition_matrix = {}
    for base in transitions:
        transition_matrix[base] = {}
        for next_base in transitions[base]:
            transition_matrix[base][next_base] = (
                transitions[base][next_base] / totals[base]
            )

    return transition_matrix


def save_to_json(matrix, filename="dna_transition_matrix.json"):
    """Save transition matrix to JSON file"""
    with open(filename, "w") as f:
        json.dump(matrix, f, indent=4)


def main():
    dna_sequence = generate_dna_sequence(50)
    print("Generated DNA sequence:")
    print(dna_sequence)

    transition_matrix = compute_transition_matrix(dna_sequence)

    print("\nTransition Matrix (Probabilities):")
    for base, transitions in transition_matrix.items():
        print(base, transitions)

    save_to_json(transition_matrix)
    print("\nTransition matrix saved to dna_transition_matrix.json")


if __name__ == "__main__":
    main()