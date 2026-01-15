import numpy as np

def predict_states(matrix, initial_vector, steps=5):
    states = [initial_vector]

    current_state = initial_vector
    for step in range(steps):
        next_state = matrix @ current_state
        states.append(next_state)
        current_state = next_state

    return states


def main():
    print("=== N-State Prediction Software ===")

    n = int(input("Enter number of states (n): "))

    print("\nEnter the transition matrix A (row by row):")
    matrix_data = []
    for i in range(n):
        row = list(map(float, input(f"Row {i+1}: ").split()))
        matrix_data.append(row)

    A = np.array(matrix_data)

    print("\nEnter the initial state vector x0:")
    x0 = np.array(list(map(float, input().split())))

    print("\n--- Prediction Results (5 Steps) ---")
    predicted_states = predict_states(A, x0, steps=5)

    for i, state in enumerate(predicted_states):
        print(f"Step {i}: {state}")


if __name__ == "__main__":
    main()
