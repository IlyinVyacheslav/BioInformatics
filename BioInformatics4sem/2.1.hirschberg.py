# Low-memory global alignment. Hirschberg’s algorithm
match = 2
mismatch = -1
gap = -2


def N_W_Score(s, t):
    score = [[0 for _ in range(len(t) + 1)] for _ in range(2)]
    for j in range(1, len(t) + 1):
        score[0][j] = gap * j
    for i in range(1, len(s) + 1):
        score[1][0] = gap * i
        for j in range(1, len(t) + 1):
            score[1][j] = max(score[0][j-1] + (match if s[i-1] == t[j-1] else mismatch), 
                              score[1][j-1] + gap, score[0][j] + gap)
        score[0] = score[1].copy()
    return score[1]


def N_W(l, s):
    ans = ["-"] * len(s)
    appended = False
    for i in range(len(s)):
        if l == s[i]:
            ans[i] = l
            appended = True
            break
    if not appended:
        ans[0] = l
    return ["".join(ans), s]            


def Hirschberg(s, t):
    n = len(s)
    m = len(t)
    if n == 0:
        return ["-" * m, t]
    elif m == 0:
        return [s, "-" * n]
    elif n == 1:
        return N_W(s, t)
    elif m == 1:
        return N_W(t, s)[::-1]
    else:
        xmid = n // 2
        left = N_W_Score(s[:xmid], t)
        right = N_W_Score(s[xmid:][::-1], t[::-1])[::-1]
        ymid = max(range(len(left)), key=lambda i : left[i] + right[i])
        opt1 = Hirschberg(s[:xmid], t[:ymid])
        opt2 = Hirschberg(s[xmid:], t[ymid:])
        return [opt1[0] + opt2[0], opt1[1] + opt2[1]]


s = input()
t = input()

res = Hirschberg(s, t)
print(res[0])
print(res[1])
