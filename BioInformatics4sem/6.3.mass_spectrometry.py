"""
Finds the peptide in the given proteome that corresponds the best with the given Spectrum
"""


def get_AA_letters_masses():
    letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    masses = [71, 103, 115, 129, 147, 57, 137, 113, 128, 113, 131, 114, 97, 128, 156, 87, 101, 99, 186, 163]
    return letters, masses


def get_XZ_letters_masses():
    letters = ['X', 'Z']
    masses = [4, 5]
    return letters, masses


def get_acid(mass, alphabet, masses):
    try:
        return alphabet[masses.index(mass)]
    except ValueError:
        return -1


def get_mass(acid, alphabet, masses):
    try:
        return masses[alphabet.index(acid)]
    except ValueError:
        return -1


def best_protein(proteome, spectral_vector, letters, masses):
    best_peptide = ""
    best_value = -1
    for i in range(len(proteome)):
        current_mass = 0
        v = 0
        j = i
        while j < len(proteome):
            current_mass += spectral_vector[v]
            v += masses[letters.index(proteome[j])]
            if v >= len(spectral_vector):
                break
            if v == len(spectral_vector) - 1 and current_mass + v + spectral_vector[v] >= best_value:
                best_value = current_mass + v + spectral_vector[v]
                best_peptide = proteome[i:j + 1]
            j += 1
    return best_peptide


def __main__():
    spectral_vector = list(map(int, input().split()))
    spectral_vector.insert(0, 0)
    proteome = input()
    if len(spectral_vector) > 30:
        letters, masses = get_AA_letters_masses()
    else:
        letters, masses = get_XZ_letters_masses()
    print(best_protein(proteome, spectral_vector, letters, masses))


if __name__ == "__main__":
    __main__()
