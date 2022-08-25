#!/bin/bash

declare -a n_arr=( 300 400 500 600 )
declare -a ant_percentage=( 2 35 50 65 80 )
declare -a alpha=( 1 2 3 4 )
declare -a beta=( 1 2 3 4 )
declare -a rho=( 0.9 0.8 0.7 0.6 )

# compile
g++ ../src/main.cpp ../src/aco.cpp -O3 -o ../src/main

echo "#ALPHA/BETA TESTS" >> "../src/results.txt"
for a in "${alpha[@]}"
do
    for b in "${beta[@]}"
    do
        echo "${a} ${b}">> "../src/results.txt"
        # generate the instance first
        python3 instance_generator.py 300 8

        # run the aco algorithm
        ../src/main ../src/input.txt $((300*8/10)) ${a} ${b} 0.8

        # compare the sequences
        python3 accuracy.py ../src/output.txt
    done
done

echo "\n#RHO TESTS" >> "../src/results.txt"
for r in "${rho[@]}"
do
    echo "${r}">> "../src/results.txt"
    # generate the instance first
    python3 instance_generator.py 300 8

    # run the aco algorithm
    ../src/main ../src/input.txt $((300*8/10)) 1 1 ${r}

    # compare the sequences
    python3 accuracy.py ../src/output.txt
done


echo "#N TEST" >> "../src/results.txt"
for n in "${n_arr[@]}"
do
    echo "${n}">> "../src/results.txt"
    # generate the instance first
    python3 instance_generator.py ${n} 8

    # run the aco algorithm
    ../src/main ../src/input.txt $((300*8/10)) ${a} ${b} 0.8

    # compare the sequences
    python3 accuracy.py ../src/output.txt
done


echo "#ANT TEST" >> "../src/results.txt"
for a in "${ant_percentage[@]}"
do
    echo "${a}">> "../src/results.txt"
    # generate the instance first
    python3 instance_generator.py 300 8

    # run the aco algorithm
    ../src/main ../src/input.txt $((300*${a}/100)) 1 1 0.8

    # compare the sequences
    python3 accuracy.py ../src/output.txt
done

