# Use the transition matrix from the JSON in order to synthesize new sequence of text based on the transition matrix

import json
import random


def load_transition_matrix(filename="word_transition_matrix.json"):
    with open(filename, "r") as f:
        data = json.load(f)
    return data["word_to_symbol"], data["transition_matrix"]


def invert_dictionary(d):
    return {v: k for k, v in d.items()}


def choose_next_state(transitions):
    states = list(transitions.keys())
    probabilities = list(transitions.values())
    return random.choices(states, probabilities)[0]


def generate_text(word_to_symbol, transition_matrix, length=50):
    symbol_to_word = invert_dictionary(word_to_symbol)

    current_symbol = random.choice(list(transition_matrix.keys()))
    generated_symbols = [current_symbol]

    for _ in range(length - 1):
        if current_symbol not in transition_matrix:
            current_symbol = random.choice(list(transition_matrix.keys()))
        else:
            current_symbol = choose_next_state(transition_matrix[current_symbol])

        generated_symbols.append(current_symbol)

    generated_words = [
        symbol_to_word[symbol] for symbol in generated_symbols
    ]

    return " ".join(generated_words)


def main():
    word_to_symbol, transition_matrix = load_transition_matrix()

    generated_text = generate_text(word_to_symbol, transition_matrix, length=50)

    print("=== Generated Text ===\n")
    print(generated_text)


if __name__ == "__main__":
    main()