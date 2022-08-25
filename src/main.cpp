#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <map>
#include <chrono>
#include "aco.h"
#include "constants.h"
#include <set>
#include <fstream>

using namespace std;
using namespace std::chrono_literals;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

void readInput(string filename, int& startSeqenceLength, int& oligoLength, string& startSequence, vector<Oligo>& spectrum) {
    ifstream infile(filename);

    int index = 0;
    infile >> startSeqenceLength;
    infile >> oligoLength;
    infile >> startSequence;

    string sequence;
    while (infile >> sequence){
        spectrum.push_back(Oligo(index, sequence));
        index++;
    }

    infile.close();
}

void writeResult(string filename, string startSequence, string finalSequence, duration<double, std::milli> main_loop_time) {
    ofstream outfile(filename);
    outfile << main_loop_time.count() << endl << startSequence << endl << finalSequence;
    outfile.close();
}

int main(int argc, char * argv[]) {
    if (argc != 6) {
        cout << "Provide all arguments <inputfile, antnumber, alpha, beta, rho>.\n";
        return 0;
    }

    int startSeqenceLength, oligoLength, antNumber = atoi(argv[2]);
    double alpha = atof(argv[3]), beta = atof(argv[4]), rho = atof(argv[5]);
    string filename = argv[1], startSequence, finalSequence;
    vector<Oligo> spectrum;
    readInput(filename, startSeqenceLength, oligoLength, startSequence, spectrum);

    ACO aco = ACO(startSeqenceLength, oligoLength, startSequence, spectrum, antNumber, alpha, beta, rho);

    cout << aco.antNumber << endl;

    auto main_loop_t1 = high_resolution_clock::now();

        aco.init();

        cout << "Init done" << endl;
        finalSequence = aco.optimize();

        cout << endl << endl;
        cout << startSequence << endl;
        cout << finalSequence << endl << endl;

    auto main_loop_t2 = high_resolution_clock::now();
    duration<double, std::milli> main_loop_time = main_loop_t2 - main_loop_t1;
    
    std::cout << "time: " << main_loop_time.count() << "ms\n";

    writeResult("../src/output.txt", startSequence, finalSequence, main_loop_time);

}


