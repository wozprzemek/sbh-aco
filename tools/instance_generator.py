import random
import copy

NUCLEOTIDES = ['A', 'C', 'G', 'T']
N = 200 # sequence length
K = 8 # slice length
NEG = int(0.06 * (N-K+1)) # negative error percentage
POS = int(0.06 * (N-K+1)) # positive error percentage

def generate_sequence():
    return [random.choice(NUCLEOTIDES) for i in range(N)]


def generate_spectrum(sequence):
    return [''.join(sequence[i:K+i]) for i in range(N-K+1)]


def write_to_file(filename, n, k, spectrum, sequence):
    with open(filename, 'w') as file:
        file.write(str(n) + '\n')
        file.write(str(k) + '\n')
        file.write(''.join(sequence) + '\n')
        for v in spectrum:
            file.write(v + '\n')


def apply_negative_error(spectrum, neg):
    # shuffle the spectrum before applying errors
    new_spectrum = copy.deepcopy(spectrum)
    random.shuffle(new_spectrum)
    # apply negative error
    if neg == 0:
        return new_spectrum
    
    print(new_spectrum[-neg:])
    return new_spectrum[:-neg]


def apply_positive_error(spectrum, pos):
    new_spectrum = copy.deepcopy(spectrum)
    appended = []
    # apply positive error
    for i in range(pos):
        v = [random.choice(NUCLEOTIDES) for i in range(K)] # single slice to be added
        new_spectrum.append(''.join(v))
        appended.append(''.join(v))
    # shuffle again
    random.shuffle(new_spectrum)
    print(appended)
    return new_spectrum

sequence = generate_sequence()
spectrum = generate_spectrum(sequence)
initial_spectrum = copy.deepcopy(spectrum)
spectrum = apply_negative_error(spectrum, NEG)
spectrum = apply_positive_error(spectrum, POS)
print("Initial spectrum:", ' '.join(initial_spectrum))
print("Initial sequence:", ''.join(sequence))
write_to_file('../src/input.txt', N, K, spectrum, sequence)