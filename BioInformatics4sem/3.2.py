weights = {
    ('A', 'U'): 1,
    ('U', 'A'): 1, 
    ('C', 'G'): 1,
    ('G', 'C'): 1,
}

s = input()
n = len(s)


def traceback(left, right, path, ans):
    if (left < right):
        cell = path[left][right]
        if cell == None:
            for k in range(left, right+1):
                ans[k] = '.'
            return
        if len(cell) == 3:
            ans[left] = '('
            ans[right] = ')'
            traceback(cell[0], cell[1], path, ans)
        elif len(cell) == 2:
            ans[left] = '.'
            ans[right] = '.'
            traceback(cell[0], cell[1], path, ans)
        else:
            traceback(cell[0], cell[1], path, ans)
            traceback(cell[2], cell[3], path, ans)
    elif left == right:
        ans[left] = '.'


dp = [[0] * n for _ in range(n)]
path = [[None] * n for _ in range(n)]


for x in range(2, n):
    for i in range(n - x):
        j = i + x
        ind = max(range(i,j), key=lambda k : dp[i][k] + dp[k+1][j])
        tmp = dp[i][ind] + dp[ind+1][j]
        diag = dp[i+1][j-1] + weights.get((s[i], s[j]), 0)
        if (diag >= tmp):
            dp[i][j] = diag
            path[i][j] = (i+1, j-1, 1) if weights.get((s[i], s[j]), 0) == 1 else (i+1,j-1)
        else:
            dp[i][j] = tmp
            path[i][j] = (i, ind, ind + 1, j)


ans = [''] * n
traceback(0, n-1, path, ans)
print(''.join(ans))