weights = {
    ('A', 'U'): 1,
    ('U', 'A'): 1, 
    ('C', 'G'): 1,
    ('G', 'C'): 1,
}

s = input()
n = len(s)

dp = [[0]* n for _ in range(n)]


for x in range(2, n):
    for i in range(n - x):
        j = i + x
        tmp = max(dp[i][k] + dp[k+1][j] for k in range(i, j))
        dp[i][j] = max(dp[i+1][j-1] + weights.get((s[i], s[j]), 0), tmp)

print(dp[0][n-1])