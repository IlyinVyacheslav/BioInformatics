import sys

sys.setrecursionlimit(10000)

n, m = map(int, input().split())
g = [[] for _ in range(n)]
for i in range(m):
    start, end = map(int, input().split())
    start -= 1
    end -= 1
    g[start].append([i, end])
visited = [False for _ in range(m)]
first = [0 for _ in range(n)]
ans = []


def euler(u):
    while first[u] < len(g[u]):
        index, v = g[u][first[u]]
        first[u] += 1
        if not visited[index]:
            visited[index] = True
            euler(v)
            ans.append(v + 1)


euler(0)
print(*ans[::-1])
