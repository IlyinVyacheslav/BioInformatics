# global alignment of 2 sequences with fine for opening the gap
from enum import Enum


class Transition(Enum):
    MAIN = 0
    TO_LOWER = 1
    TO_UPPER = 2
    TO_MAIN = 3
    RIGHT = 4
    DOWN = 5


def generate_list(x, y, value=0, transition=None):
    return [[[value, transition] for _ in range(x)] for _ in range(y)]


mismatch = -1
match = 1
INF = -10 ** 6

s = input()
t = input()
gap_start, gap_continue = map(int, input().split())

s_len, t_len = len(s) + 1, len(t) + 1
main = generate_list(s_len, t_len)
upper = generate_list(s_len, t_len, value=0, transition=Transition.DOWN)
lower = generate_list(s_len, t_len, value=0, transition=Transition.RIGHT)
upper[0][0], lower[0][0] = [INF, Transition.DOWN], [INF, Transition.RIGHT]
main[0][0] = [0, Transition.MAIN]
for i in range(1, s_len):
    gap_penalty = gap_start + gap_continue * i
    lower[0][i] = [gap_penalty, Transition.RIGHT]
    upper[0][i] = [INF, Transition.DOWN]
    main[0][i] = [gap_penalty, Transition.RIGHT]

for i in range(1, t_len):
    gap_penalty = gap_start + gap_continue * i
    lower[i][0] = [INF, Transition.RIGHT]
    upper[i][0] = [gap_penalty, Transition.DOWN]
    main[i][0] = [gap_penalty, Transition.DOWN]

for i in range(1, t_len):
    for j in range(1, s_len):
        upper[i][j] = [upper[i - 1][j][0] + gap_continue, Transition.DOWN] if upper[i - 1][j][0] + gap_continue >= \
                                                                              main[i - 1][j][
                                                                                  0] + gap_start + gap_continue else [
            main[i - 1][j][0] + gap_start + gap_continue, Transition.TO_MAIN]
        lower[i][j] = [lower[i][j - 1][0] + gap_continue, Transition.RIGHT] if lower[i][j - 1][0] + gap_continue >= \
                                                                               main[i][j - 1][
                                                                                   0] + gap_start + gap_continue else [
            main[i][j - 1][0] + gap_start + gap_continue, Transition.TO_MAIN]
        main_cell_value = main[i - 1][j - 1][0] + (match if s[j - 1] == t[i - 1] else mismatch)
        main_cell = max(main_cell_value, upper[i][j][0], lower[i][j][0])
        main[i][j] = [main_cell, Transition.MAIN if main_cell == main_cell_value else (
            Transition.TO_UPPER if main_cell == upper[i][j][0] else Transition.TO_LOWER)]

i, j = t_len - 1, s_len - 1
s_aligned, t_aligned = [], []

flag = Transition.MAIN

while i > 0 or j > 0:
    current_cell = main[i][j]
    if flag == Transition.DOWN:
        current_cell = upper[i][j]
    elif flag == Transition.RIGHT:
        current_cell = lower[i][j]
    transition = current_cell[1]

    if transition == Transition.MAIN and flag == Transition.MAIN:
        i -= 1
        j -= 1
        s_aligned.append(s[j])
        t_aligned.append(t[i])
    elif transition == Transition.TO_UPPER:
        flag = Transition.DOWN
    elif transition == Transition.DOWN:
        i -= 1
        s_aligned.append('-')
        t_aligned.append(t[i])
    elif transition == Transition.TO_LOWER:
        flag = Transition.RIGHT
    elif transition == Transition.RIGHT:
        j -= 1
        s_aligned.append(s[j])
        t_aligned.append('-')
    elif transition == Transition.TO_MAIN and flag == Transition.DOWN:
        i -= 1
        s_aligned.append('-')
        t_aligned.append(t[i])
        flag = Transition.MAIN
    elif transition == Transition.TO_MAIN and flag == Transition.RIGHT:
        j -= 1
        s_aligned.append(s[j])
        t_aligned.append('-')
        flag = Transition.MAIN
    else:
        break

print(''.join(reversed(s_aligned)))
print(''.join(reversed(t_aligned)))
