#!/bin/bash

# generate the instance first
python3 instance_generator.py

# run the aco algorithm
g++ ../src/main.cpp -O3 -o ../src/main
../src/main ../src/input.txt

# plot the results
# python3 plot_results.py