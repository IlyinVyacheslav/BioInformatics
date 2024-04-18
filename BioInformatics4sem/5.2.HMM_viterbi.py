"""
Returns the most probable sequence of hidden states for the given nucleotide sequence
based on Hidden Markov Model using Viterbi algorithm
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
dp = [[0.0 for _ in range(len(s))] for _ in range(2)]
dp[0][0] = 0.5 * hmm.b(Hidden.P, getattr(Observable, s[0]))
dp[1][0] = 0.5 * hmm.b(Hidden.N_P, getattr(Observable, s[0]))
ans = [["" for _ in range(len(s))] for _ in range(2)]
for t in range(1, len(s)):
    dp[0][t], ans[0][t] = get_max(
        dp[0][t - 1] * hmm.a(Hidden.P, Hidden.P) * hmm.b(Hidden.P, getattr(Observable, s[t])),
        dp[1][t - 1] * hmm.a(Hidden.N_P, Hidden.P) * hmm.b(Hidden.P, getattr(Observable, s[t])))
    dp[1][t], ans[1][t] = get_max(
        dp[0][t - 1] * hmm.a(Hidden.P, Hidden.N_P) * hmm.b(Hidden.N_P, getattr(Observable, s[t])),
        dp[1][t - 1] * hmm.a(Hidden.N_P, Hidden.N_P) * hmm.b(Hidden.N_P, getattr(Observable, s[t])))

hidden_states = ["" for _ in range(len(s))]
hidden_states[len(s) - 1], pos = ("P", 0) if dp[0][len(s) - 1] >= dp[1][len(s) - 1] else ("N", 1)
for i in range(len(s) - 2, -1, -1):
    hidden_states[i] = ans[pos][i + 1]
    pos = 0 if hidden_states[i] == "P" else 1
print(*hidden_states, sep="")
