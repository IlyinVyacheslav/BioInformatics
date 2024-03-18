# returns the weight of global alignment of 2 sequences with fine for opening the gap
from enum import Enum


class Shift(Enum):
    default = -1
    diag_match = 0
    diag_skip = 1
    right = 2
    down = 3


def generate_list(x, y):
    return [[0 for i in range(x)] for j in range(y)]


mismatch = -1
match = 1
INF = -10 ** 6

s = input()
t = input()
gap_start, gap_continue = map(int, input().split())

s_len, t_len = len(s) + 1, len(t) + 1
main = generate_list(s_len, t_len)
upper = generate_list(s_len, t_len)
lower = generate_list(s_len, t_len)
upper[0][0], lower[0][0] = INF, INF
for i in range(1, s_len):
    gap_penalty = gap_start + gap_continue * i
    lower[0][i] = gap_penalty
    upper[0][i] = INF
    main[0][i] = gap_penalty

for i in range(1, t_len):
    gap_penalty = gap_start + gap_continue * i
    lower[i][0] = INF
    upper[i][0] = gap_penalty
    main[i][0] = gap_penalty

for i in range(1, t_len):
    for j in range(1, s_len):
        upper[i][j] = max(upper[i - 1][j] + gap_continue, main[i - 1][j] + gap_start + gap_continue)
        lower[i][j] = max(lower[i][j - 1] + gap_continue, main[i][j - 1] + gap_start + gap_continue)
        main[i][j] = max(main[i - 1][j - 1] + (match if s[j - 1] == t[i - 1] else mismatch),
                         upper[i][j], lower[i][j])

print(main[t_len - 1][s_len - 1])
