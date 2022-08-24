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

void writeResult(string filename, string startSequence, string finalSequence) {
    ofstream outfile(filename);
    outfile << startSequence << endl << finalSequence;
    outfile.close();
}

int main(int argc, char * argv[]) {

    if (argc != 4) {
        cout << "Provide all arguments <inputfile, iterations, antnumber>.\n";
        return 0;
    }

    int startSeqenceLength, oligoLength, iterations = atoi(argv[2]), antNumber = atoi(argv[3]);
    string filename = argv[1], startSequence, finalSequence;
    vector<Oligo> spectrum;
    readInput(filename, startSeqenceLength, oligoLength, startSequence, spectrum);

    ACO aco = ACO(startSeqenceLength, oligoLength, startSequence, spectrum, iterations, antNumber);

    cout << aco.antNumber << endl;

    aco.init();

    cout << "Init done" << endl;
    finalSequence = aco.optimize();

    cout << endl << endl;
    cout << startSequence << endl;
    cout << finalSequence << endl << endl;

    writeResult("../src/output.txt", startSequence, finalSequence);

}


