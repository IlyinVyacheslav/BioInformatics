"""
Calculates theoretical ideal spectrum of peptide
"""
alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
masses = [71, 103, 115, 129, 147, 57, 137, 113, 128, 113, 131, 114, 97, 128, 156, 87, 101, 99, 186, 163]


def get_mass(letter):
    return masses[alphabet.index(letter)]


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


print(*get_ideal_spectrum(input()))
