import random
from collections import defaultdict, deque

S = "ACACCAAATGGAACGAAAATCAGAACCCTAGAATGTTTTTGGCCATGATCACATATATAACCAGAAATCAGCCCGAATGGTTCAGAAATATTCTAAGTATTGCTCCAATAATGTTCTCAAACAAAATGGCGAGACTAGGTAAGGGGTACATGTTTGAAAGCAAGAGTATGAAACTTAGAACTCAAATACCTGCAGAAATGCTAGCAAACATCGATTTGAGGTATTTCAATGATTCAACAAGAAAGAAAATTGAAAAAATCCGGCCACTCTTAATAGATGGAACTGCATCACTGAGCCCTGGAATGATGATGGGCATGTTCAATATGTTAAGCACTGTATTGGGCGTCTCCATCCTGAATCTTGGACAAAAGAGATACACCAAGACTACTTACTGGTGGGATGGTCTTCAATCGTCTGATGATTTTGCTCTGATTGTGAATGCACCCAACTATGCAGGAATTCAAGCTGGAGTTGACAGGTTTTATCGAACCTGTAAGCTGCTCGGAATCAACATGAGCAAGAAAAAGTCCTACATAAACAGAACAGGTACTTTTGAATTCACGAGTTTTTTCTATCGTTATGGGTTTGTTGCCAATTTCAGCATGGAGCTTCCCAGTTTTGGGGTGTCTGGGATCAACGAGTCTGCAGACATGAGTATTGGAGTCACTGTCATCAAAAACAATATGATAAACAATGATCTTGGTCCAGCAACCGCTCAAATGGCCCTTCAGTTATTCATCAAAGATTACAGGTACACATATCGATGCCACAGAGGTGACACACAAATACAAACCCGAAGATCATTTGAAATAAAGAAACTGTGGGATCAAACCCGTTCCAAAACTGGGCTGCTGGTCTCTGATGGAGGTCCTAATTTGTACAACATTAGAAATCTCCACATTCCTGAAGTCTGCTTGAAATGGGATTTGATGGATGAGGATTACCAGGGGCGTTTATGCAACCCATTGAACCCATTTGTCAGTCATAAAGAGATTGAATCAGTGAACAATGCAGTGATGATGCCGGCACATGGTCCAGCCAAAATTATGGAGTATGATGCTGTTGCAACAACACACTCCTGGGTCCCCAAAAGGAATCGATCCATCTTGAATACGAGCCAAAGGGGAATACTTGAGGATGAACAAATGTATCAGAGGTGCTGCAATTTATTTGAAAAATTCTTCCCCAGCAGTTCATACAGAAGACCAGTCGGAATATCCAGTATGGTAGAGGCTATGATTTCCAGAGCCCGAATTGATGCACGGATTGATTTCGAATCTGGAAGGATAAAGAAAGAGGAGTTCACTGAGATCATGAAGATCTGTTCCACCATTGAAGAACTCAGACGGCAAAAATAGGGAATTTGGCTTGTCCTTCGTGAAAAAATGCCTTGTTTCTACT"

NR_SAMPLES = 2000
SAMPLE_LENGTH = 100
K = 90

reads = [S[random.randint(0, len(S) - SAMPLE_LENGTH):][:SAMPLE_LENGTH] for _ in range(NR_SAMPLES)]

def build_debruijn(reads, k):
    graph = defaultdict(list)
    indeg, outdeg = defaultdict(int), defaultdict(int)
    for read in reads:
        for i in range(len(read) - k + 1):
            prefix = read[i:i + k - 1]
            suffix = read[i + 1:i + k]
            graph[prefix].append(suffix)
            outdeg[prefix] += 1
            indeg[suffix] += 1
    return graph, indeg, outdeg

def eulerian_path(graph, indeg, outdeg):
    start = None
    for node in graph:
        if outdeg[node] - indeg[node] == 1:
            start = node
            break
    if not start:
        start = next(iter(graph))

    stack, path = [start], []
    while stack:
        v = stack[-1]
        if graph[v]:
            u = graph[v].pop()
            stack.append(u)
        else:
            path.append(stack.pop())
    path.reverse()
    return path

def assemble(reads, k):
    graph, indeg, outdeg = build_debruijn(reads, k)
    path = eulerian_path(graph, indeg, outdeg)
    assembled = path[0] + ''.join(node[-1] for node in path[1:])
    return assembled

assembled_seq = assemble(reads, K)

print(f"Original length: {len(S)}")
print(f"Reconstructed length: {len(assembled_seq)}")

if S in assembled_seq or assembled_seq in S:
    print("Successfully reconstructed the sequence (or close match)")
else:
    print("Partial reconstruction — try higher NR_SAMPLES or smaller K")

print("Reconstructed start:", assembled_seq[:200])
