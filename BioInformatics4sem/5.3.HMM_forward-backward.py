"""
Returns the probabilities that position is in promoter or non-promoter state
for the given nucleotide sequence based on Hidden Markov Model using Forward-Backward algorithm
"""
from enum import Enum, auto


class Hidden(Enum):
    P = auto()
    N_P = auto()


class Observable(Enum):
    C = auto()
    G = auto()
    A = auto()
    T = auto()


class HMM:
    def __init__(self):
        self.probability_matrix = [
            # C,   G,   A,   T
            [0.4, 0.4, 0.1, 0.1],
            [0.2, 0.2, 0.3, 0.3]
        ]
        self.hidden_state_matrix = [
            # P    NP
            [0.9, 0.1],
            [0.2, 0.8]
        ]

    def b(self, hidden_state, observable_state):
        return self.probability_matrix[hidden_state.value - 1][observable_state.value - 1]

    def a(self, hs_1, hs_2):
        return self.hidden_state_matrix[hs_1.value - 1][hs_2.value - 1]


def get_max(a, b):
    if a >= b:
        return a, "P"
    else:
        return b, "N"


hmm = HMM()
s = input()
A = [[0.0 for _ in range(len(s))] for _ in range(2)]
B = [[0.0 for _ in range(len(s))] for _ in range(2)]
P = [[0.0 for _ in range(len(s))] for _ in range(2)]
A[0][0] = 0.5 * hmm.b(Hidden.P, getattr(Observable, s[0]))
A[1][0] = 0.5 * hmm.b(Hidden.N_P, getattr(Observable, s[0]))

for t in range(1, len(s)):
    A[0][t] = (A[0][t - 1] * hmm.a(Hidden.P, Hidden.P) * hmm.b(Hidden.P, getattr(Observable, s[t])) +
               A[1][t - 1] * hmm.a(Hidden.N_P, Hidden.P) * hmm.b(Hidden.P, getattr(Observable, s[t])))
    A[1][t] = (A[0][t - 1] * hmm.a(Hidden.P, Hidden.N_P) * hmm.b(Hidden.N_P, getattr(Observable, s[t])) +
               A[1][t - 1] * hmm.a(Hidden.N_P, Hidden.N_P) * hmm.b(Hidden.N_P, getattr(Observable, s[t])))
B[0][len(s) - 1], B[1][len(s) - 1] = 1, 1
for t in range(len(s) - 2, -1, -1):
    B[0][t] = (B[0][t + 1] * hmm.a(Hidden.P, Hidden.P) * hmm.b(Hidden.P, getattr(Observable, s[t + 1])) +
               B[1][t + 1] * hmm.a(Hidden.P, Hidden.N_P) * hmm.b(Hidden.N_P, getattr(Observable, s[t + 1])))
    B[1][t] = (B[0][t + 1] * hmm.a(Hidden.N_P, Hidden.P) * hmm.b(Hidden.P, getattr(Observable, s[t + 1])) +
               B[1][t + 1] * hmm.a(Hidden.N_P, Hidden.N_P) * hmm.b(Hidden.N_P, getattr(Observable, s[t + 1])))

P_sum = A[0][len(s) - 1] + A[1][len(s) - 1]
for t in range(len(s)):
    P[0][t] = A[0][t] * B[0][t] / P_sum
    P[1][t] = A[1][t] * B[1][t] / P_sum

for row in P:
    print(' '.join('{:.2f}'.format(x) for x in row))
