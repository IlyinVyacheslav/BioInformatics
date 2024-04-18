"""
Finds the peptide that could have generated given ideal Spectrum
"""
alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
masses = [71, 103, 115, 129, 147, 57, 137, 113, 128, 113, 131, 114, 97, 128, 156, 87, 101, 99, 186, 163]


def get_mass(letter):
    return masses[alphabet.index(letter)]


def get_acid(mass):
    try:
        return alphabet[masses.index(mass)]
    except ValueError:
        return -1


def build_graph(Spectrum):
    n = len(Spectrum)
    edges = [['' for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            acid = get_acid(Spectrum[j] - Spectrum[i])
            if acid != -1:
                edges[i][j] = acid
    return edges


def dfs(graph, v, end, path, all_paths):
    path.append(v)
    if v == end:
        all_paths.append(path.copy())
    else:
        for neighbour, acid in enumerate(graph[v]):
            if acid and neighbour not in path:
                dfs(graph, neighbour, end, path, all_paths)
    path.pop()


def get_ideal_spectrum(peptide):
    res = []
    curr = 0
    for c in peptide:
        res.append(curr)
        curr += get_mass(c)

    curr = 0
    for c in peptide[::-1]:
        curr += get_mass(c)
        res.append(curr)
    return sorted(set(res))


def decoding_ideal_spectrum(spectrum):
    graph = build_graph(spectrum)
    all_paths = []
    dfs(graph, 0, len(spectrum) - 1, [], all_paths)
    for path in all_paths:
        peptide = "".join([graph[path[i]][path[i+1]] for i in range(len(path) - 1)])
        if get_ideal_spectrum(peptide).__eq__(spectrum):
            return peptide


def __main__():
    Spectrum = list(map(int, input().split()))
    print(decoding_ideal_spectrum(Spectrum))


if __name__ == "__main__":
    __main__()
