from operator import le
import random
import copy

NUCLEOTIDES = ['A', 'C', 'G', 'T']
N = 200 # sequence length
K = 8 # slice length
NEG = int(0.06 * (N-K+1)) # negative error percentage
POS = int(0.1 * (N-K+1)) # positive error percentage

def generate_sequence():
    return [random.choice(NUCLEOTIDES) for i in range(N)]


def generate_spectrum(sequence):
    return [''.join(sequence[i:K+i]) for i in range(N-K+1)]

def write_to_file(filename, n, k, spectrum):
    with open(filename, 'w') as file:
        file.write(str(n) + '\n')
        file.write(str(k) + '\n')
        for v in spectrum:
            file.write(v + '\n')


def apply_negative_error(spectrum, neg):
    # shuffle the spectrum before applying errors
    new_spectrum = copy.deepcopy(spectrum)
    random.shuffle(new_spectrum)
    # apply negative error
    return new_spectrum[:-neg]


def apply_positive_error(spectrum, pos):
    new_spectrum = copy.deepcopy(spectrum)
    # apply positive error
    for i in range(pos):
        v = [random.choice(NUCLEOTIDES) for i in range(K)] # single slice to be added
        new_spectrum.append(''.join(v))
    # shuffle again
    random.shuffle(new_spectrum)
    return new_spectrum

sequence = generate_sequence()
spectrum = generate_spectrum(sequence)

apply_negative_error(spectrum, NEG)
apply_positive_error(spectrum, POS)

write_to_file('input.txt', N, K, spectrum)