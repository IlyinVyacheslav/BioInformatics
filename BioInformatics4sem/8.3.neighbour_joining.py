# returns phylogenetic tree in Newick format from the given distance matrix using UPGMA
import sys

class Cluster(object):
    def __init__(self, name):
        self.name = name

    def merge(self, other, d1, d2):
        self.name = f"({self.name}:{d1:.2f},{other.name}:{d2:.2f})"


class DistanceTable(object):
    def __init__(self, D):
        self.D = D
        self.clusters = [Cluster(chr(ord('A') + i)) for i in range(len(D))]

    def find_nearest_vertices(self):
        min_pair = None
        N = len(self.D)
        min_Q = float('inf')
        for i in range(N):
            for j in range(i + 1, N):
                Q_ij = (N - 2) * self.D[i][j] - sum(self.D[i][k] + self.D[j][k] for k in range(N))
                if Q_ij < min_Q:
                    min_Q = Q_ij
                    min_pair = (i, j)
        return min_pair

    def is_reducible(self):
        return len(self.D) >= 2

    def merge_clusters(self, indices):
        i, j = indices
        N = len(self.D)
        if N > 2:
            sum_i = sum(self.D[i][k] for k in range(N))
            sum_j = sum(self.D[j][k] for k in range(N))
            delta_iu = (self.D[i][j] + (sum_i - sum_j) / (N - 2)) / 2
            delta_ju = self.D[i][j] - delta_iu
        else:
            delta_iu = delta_ju = self.D[i][j] / 2

        for k in range(N):
            if k != i and k != j:
                self.D[k][i] = self.D[i][k] = (self.D[i][k] + self.D[j][k] - self.D[i][j]) / 2

        self.clusters[i].merge(self.clusters[j], delta_iu, delta_ju)
        self.clusters.pop(j)

        del self.D[j]
        for k in range(len(self.D)):
            del self.D[k][j]

    def print(self):
        max_width = max(len(str(cluster.name)) for cluster in self.clusters)
        max_width = max(max_width, 5)

        print(" " * (max_width + 3), end="")
        for cluster in self.clusters:
            print(f"{str(cluster.name).ljust(max_width)} ", end="")
        print()

        for i, row in enumerate(self.D):
            print(f"{str(self.clusters[i].name).ljust(max_width)} ", end="")
            for value in row:
                print(f"{str(value).ljust(max_width)} ", end="")
            print()


def build_distance_matrix(n, sequences):
    D = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            D[i][j] = D[j][i] = float(sum(c1 != c2 for c1, c2 in zip(sequences[i], sequences[j])))
    return D


def __main__():
    input_List = input().split()
    n = int(input_List[0])
    SEQ = input_List[1:]
    D = build_distance_matrix(n, SEQ)
    Table = DistanceTable(D)
    while Table.is_reducible():
        Table.merge_clusters(Table.find_nearest_vertices())
    print(Table.clusters[0].name)


if __name__ == "__main__":
    __main__()
