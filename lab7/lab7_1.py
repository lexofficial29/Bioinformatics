# Download from NCBI a seq of your choosing, design an application which is able to detect repetitions of patterns with sizes between 3b and 6b. 
# Note that this repetitions must be found one after the other.
# ex: AAT => AATAATAAT...

def find_repeats_everywhere(seq, pattern):
    seq = seq.upper()
    pattern = pattern.upper()

    size = len(pattern)
    results = []

    for i in range(len(seq) - size):

        repeats = 0
        pos = i

        while seq[pos:pos + size] == pattern:
            repeats += 1
            pos += size

        if repeats > 2:
            results.append({
                "start_index": i,
                "repeat_count": repeats,
                "full_repeat_sequence": seq[i:i + repeats * size]
            })

    return results

seq = "ATGCTAACGGTTCGATACCGGATTTGACCTTAGGCTTACGGAATAATAATGCTTACCGATAGGCTTAACGGTAGCTTACCGATGGAATTCGATGCCGTTAACCGGTAGCTAACGGTTCGATAGGCGTAACCTTAGGCTAACGGAATGCTTACCGATAGGCTTACGGTAGCTAACCGATGGAATTCGATGGCTAACCGAATGCTTAGGCTAACGGTAGCTTACCGATAGGCTTACCGGAATGCTTACCGATACCGGA"

pattern = input("Enter pattern (3–6 bases): ").strip()

if not (3 <= len(pattern) <= 6):
    print("Pattern must be 3 to 6 bases long.")
else:
    repetitions = find_repeats_everywhere(seq, pattern)

for r in repetitions:
    print(f"\nFound at index {r['start_index']}:")
    print(f"Repeats: {r['repeat_count']}")
    print(f"Full sequence: {r['full_repeat_sequence']}")

