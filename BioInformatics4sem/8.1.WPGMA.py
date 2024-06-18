# returns phylogenetic tree in Newick format from the given distance matrix using WPGMA

class Cluster(object):
    def __init__(self, name):
        self.name = name
        self.height = 0

    def merge(self, other, distance):
        self.name = f"({self.name}:{distance - self.height:.2f},{other.name}:{distance - other.height:.2f})"
        self.height = distance


class DistanceTable(object):
    def __init__(self, D):
        self.distances = D
        self.clusters = [Cluster(chr(ord('A') + i)) for i in range(len(D))]

    def find_minimum(self):
        pair = None
        min_v = float('inf')
        for i in range(len(self.distances)):
            for j in range(i + 1, len(self.distances)):
                if self.distances[i][j] < min_v:
                    min_v = self.distances[i][j]
                    pair = (i, j)
        return pair

    def is_reducible(self):
        return len(self.distances) >= 2

    def merge_clusters(self, indices):
        first, second = indices
        self.clusters[first].merge(self.clusters[second], self.distances[first][second] / 2)
        self.clusters.pop(second)

        for i in range(len(self.distances)):
            self.distances[i][first] = self.distances[first][i] = (self.distances[i][first] + self.distances[i][
                second]) / 2
        del self.distances[second]
        for i in range(len(self.distances)):
            del self.distances[i][second]

    def print(self):
        max_width = max(len(str(cluster.name)) for cluster in self.clusters)
        max_width = max(max_width, 5)

        print(" " * (max_width + 3), end="")
        for cluster in self.clusters:
            print(f"{str(cluster.name).ljust(max_width)} ", end="")
        print()

        for i, row in enumerate(self.distances):
            print(f"{str(self.clusters[i].name).ljust(max_width)} ", end="")
            for value in row:
                print(f"{str(value).ljust(max_width)} ", end="")
            print()


def __main__():
    n = int(input())
    distances = [list(map(float, input().split())) for _ in range(n)]
    Table = DistanceTable(distances)
    while Table.is_reducible():
        Table.merge_clusters(Table.find_minimum())
    print(Table.clusters[0].name)


if __name__ == "__main__":
    __main__()
