mu = -1
sigma = -2
match = 1

s = input()
t = input()
s_len = len(s) + 1
t_len = len(t) + 1
dp = [[0 for i in range(s_len)] for j in range(t_len)]
for i in range(s_len):
    dp[0][i] = sigma * i
for i in range(t_len):
    dp[i][0] = sigma * i

for i in range(1, t_len):
    for j in range(1, s_len):
        dp[i][j] = max(dp[i - 1][j - 1] + match if s[j - 1] == t[i - 1] else dp[i - 1][j - 1] + mu,
                       dp[i - 1][j] + sigma, dp[i][j - 1] + sigma)

i = t_len - 1
j = s_len - 1
s_aligned = str()
t_aligned = str()

while True:
    if i == 0 and j == 0:
        break
    elif i >= 1 and dp[i - 1][j] + sigma == dp[i][j]:
        s_aligned += "-"
        t_aligned += t[i - 1]
        i -= 1
    elif j >= 1 and dp[i][j - 1] + sigma == dp[i][j]:
        s_aligned += s[j - 1]
        t_aligned += "-"
        j -= 1
    elif i >= 1 and j >= 1 and (dp[i - 1][j - 1] + match == dp[i][j] or dp[i - 1][j - 1] + mu == dp[i][j]):
        s_aligned += s[j - 1]
        t_aligned += t[i - 1]
        i -= 1
        j -= 1

print(s_aligned[::-1])
print(t_aligned[::-1])
