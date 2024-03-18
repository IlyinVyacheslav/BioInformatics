# global alignment of 2 sequences using BLOSUM62 matrix
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")

sigma = -5

s = input()
t = input()
s_len, t_len = len(s) + 1, len(t) + 1
dp = [[0 for i in range(s_len)] for j in range(t_len)]
for i in range(s_len):
    dp[0][i] = sigma * i
for i in range(t_len):
    dp[i][0] = sigma * i

for i in range(1, t_len):
    for j in range(1, s_len):
        dp[i][j] = max(dp[i - 1][j - 1] + blosum62[s[j-1]][t[i-1]],
                       dp[i - 1][j] + sigma, dp[i][j - 1] + sigma)

i, j = t_len - 1, s_len - 1
s_aligned, t_aligned = str(), str()

while True:
    if i == 0 and j == 0:
        break
    elif i >= 1 and dp[i - 1][j] + sigma == dp[i][j]:
        i -= 1
        s_aligned += "-"
        t_aligned += t[i]
    elif j >= 1 and dp[i][j - 1] + sigma == dp[i][j]:
        j -= 1
        s_aligned += s[j]
        t_aligned += "-"
    else:
        i -= 1
        j -= 1
        s_aligned += s[j]
        t_aligned += t[i]

print(s_aligned[::-1])
print(t_aligned[::-1])
