#!/bin/bash

n_arr=( 300 400 500 600 700 )
ant_percentage=( 20 35 50 65 80 )
repeats=3 # repeats per one measurement

for i in {1..5}
do
    # generate the instance first
    python3 instance_generator.py 300 8

    # run the aco algorithm
    g++ ../src/main.cpp ../src/aco.cpp -O3 -o ../src/main
    ../src/main ../src/input.txt 100 100

    # compare the sequences
    python3 accuracy.py ../src/output.txt
done