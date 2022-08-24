#!/bin/bash

for i in {1..10}
do
    # generate the instance first
    python3 instance_generator.py

    # run the aco algorithm
    # g++ ../src/main.cpp -O3 -o ../src/main
    g++ ../src/main.cpp ../src/aco.cpp -O3 -o ../src/main
    ../src/main ../src/input.txt 100 100

    # compare the sequences
    python3 accuracy.py ../src/output.txt
done