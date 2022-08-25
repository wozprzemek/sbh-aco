import Levenshtein as lev
import sys

filename = sys.argv[1]

data = []

with open(filename) as infile:
    data = infile.read().splitlines()

accuracy = lev.ratio(data[1], data[2])

with open('../src/results.txt', 'a') as outfile:
    outfile.write(str(accuracy) + " " + data[0] + " " + data[1] + " " + data[2] + "\n")


