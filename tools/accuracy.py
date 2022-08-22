import Levenshtein as lev
import sys

filename = sys.argv[1]

sequences = []

with open(filename) as infile:
    sequences = infile.read().splitlines()

accuracy = lev.ratio(sequences[0], sequences[1])

with open('../src/results.txt', 'a') as outfile:
    outfile.write(str(accuracy) + " " + sequences[0] + " " + sequences[1] + "\n")


