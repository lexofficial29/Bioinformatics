# 2. Use a random english text of arout 300 characters (that implies spaces and punctuation too) and compute the transition probabilities between words.
# Store the transition matrix as a JSON file. For ease of implamantation you could represent each unique word by using a symbol of your choice (ASCII).

import random
import json
import string
from collections import defaultdict


def generate_random_text(length=300):
    words = [
        "the", "cat", "sat", "on", "mat", "and", "dog", "ran", "away",
        "because", "it", "was", "happy", "sad", "big", "small", "quick",
        "slow", "blue", "green", "tree", "house", "street", "car", "music"
    ]

    punctuation = [".", ",", "!", "?"]

    text = ""
    while len(text) < length:
        word = random.choice(words)
        text += word
        if random.random() < 0.15:
            text += random.choice(punctuation)
        text += " "

    return text[:length]


def tokenize_text(text):
    translator = str.maketrans("", "", string.punctuation)
    clean_text = text.lower().translate(translator)
    return clean_text.split()


def map_words_to_symbols(words):
    unique_words = sorted(set(words))
    symbol_map = {}

    for i, word in enumerate(unique_words):
        symbol_map[word] = chr(65 + i)

    return symbol_map


def compute_transition_matrix(words, symbol_map):
    transitions = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)

    for i in range(len(words) - 1):
        current_symbol = symbol_map[words[i]]
        next_symbol = symbol_map[words[i + 1]]
        transitions[current_symbol][next_symbol] += 1
        totals[current_symbol] += 1

    matrix = {}
    for symbol in transitions:
        matrix[symbol] = {}
        for next_symbol in transitions[symbol]:
            matrix[symbol][next_symbol] = (
                transitions[symbol][next_symbol] / totals[symbol]
            )

    return matrix


def main():
    text = generate_random_text(300)
    print("Generated Text:\n")
    print(text)

    words = tokenize_text(text)
    symbol_map = map_words_to_symbols(words)

    transition_matrix = compute_transition_matrix(words, symbol_map)

    output = {
        "word_to_symbol": symbol_map,
        "transition_matrix": transition_matrix
    }

    with open("word_transition_matrix.json", "w") as f:
        json.dump(output, f, indent=4)

    print("\nTransition matrix saved to word_transition_matrix.json")


if __name__ == "__main__":
    main()