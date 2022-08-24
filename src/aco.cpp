#include <vector>
#include <set>
#include <cstdlib>
#include "constants.h"
#include "aco.h"
#include <string>
#include <math.h>
#include <iostream>
#include <chrono>
#include <time.h>

using namespace std::chrono_literals;

Ant::Ant(int oligoCount, int startingIndex) {
    this->currentOligoIndex = startingIndex; //setting the starting point
    this->path.push_back(startingIndex);

    //filling the notVisited set
    for (int i = 0; i < oligoCount; i++) {
        if (i != startingIndex) {
            notVisited.insert(i);
        }
    }
}

void Ant::reset(int oligoCount, int startingIndex) {
    this->currentOligoIndex = startingIndex;
    this->path.clear();
    this->path.push_back(currentOligoIndex);
    this->notVisited.clear();
    for (int i = 0; i < oligoCount; i++) {
        if (i != startingIndex) {
            notVisited.insert(i);
        }
    }

}

void Ant::move(std::vector<std::vector<double>> T, std::vector<std::vector<double>> oligoVisibilities) {
    double random = ((double)rand() / (RAND_MAX));
    double probability;
    double pN; //numerator probability
    double pD = 0.0; //denominator probability
    double sum = 0;

    for (std::set<int>::iterator j = notVisited.begin(); j != notVisited.end(); j++) {
        pD += (pow(T[currentOligoIndex][*j], ALPHA)) * (pow(oligoVisibilities[currentOligoIndex][*j], BETA)); //denominator probability
    }

    for (std::set<int>::iterator i = notVisited.begin(); i != notVisited.end(); i++) { //iterate over not visited cities
        /*PART 1--------calculate probability for a oligo-----------*/
        pN = (pow(T[currentOligoIndex][*i], ALPHA)) * (pow(oligoVisibilities[currentOligoIndex][*i], BETA)); //numerator probability

        probability = pN / pD;

        /*PART 2--------try to select the oligo for which the probability was calculated-----------*/
        random -= probability;
        if (random <= 0) {
            // cout << "random after subtraction " << random << endl;
            currentOligoIndex = *i;
            path.push_back(currentOligoIndex); //add selected oligo to ant's path
            notVisited.erase(i); //remove selected oligo from not visited
            break;
        }
    }
}

Oligo::Oligo(int index, std::string sequence) {
    this->index = index;
    this->sequence = sequence;
}

ACO::ACO(int startSeqenceLength, int oligoLength, std::string startSequence, std::vector<Oligo> spectrum, int iterations, int antNumber) {
    this->startSeqenceLength = startSeqenceLength;
    this->oligoLength = oligoLength;
    this->startSequence = startSequence;
    this->spectrum = spectrum;
    this->iterations = iterations;
    this->antNumber = antNumber;
}

void ACO::init() {
    this->oligoCount = spectrum.size();
    this->startingOligo = this->startSequence.substr(0, this->oligoLength);
    this->startingOligoIndex = this->findStartingOligo(this->startingOligo);
    this->minimum = this->oligoCount * this->oligoLength;
    this->calculateOverlaps(); //initialize oligo overlap table
    this->calculateVisibilities(); //initialize oligo visibilities table
    this->initializePheromones(); //initialize pheromones table

    this->spawnAnts(); //spawn all ants
}

int ACO::overlap(std::string a, std::string b, int len) {
    int overlapScore = len;
    for (int i=1; i<len+1; i++) {
        std::string aEnd =  a.substr(len - i, i);
        std::string bStart = b.substr(0, i);

        if (aEnd == bStart) {
            overlapScore = len - i; // lower overlapScore == bigger overlap
        }
    }
    
    return overlapScore;
}

double ACO::calculateLk(Ant ant) {
    float Lk = 0;
    for (int i = 1; i < this->oligoCount; i++) {
        Lk += this->oligoOverlaps[ant.path[i - 1]][ant.path[i]];
    }
    return Lk;
}

void ACO::updateTrailLevels() {
    for (int l = 0; l < this->antNumber; l++) { //iterate over ants
        double Lk = this->calculateLk(ants[l]);
        for (int i = 1; i < this->oligoCount; i++) { //add the pheromones along the ant's path
            this->T[this->ants[l].path[i - 1]][this->ants[l].path[i]] += Q / Lk; 
        }
        if (Lk < this->minimum) {
            this->minimum = Lk; //check if the path is the shortest yet
            this->shortestPath = this->ants[l].path;
        }
    }
}


void ACO::evaporateTrailLevels() {
    for (int i = 0; i < oligoCount; i++) {
        for (int j = 0; j < oligoCount; j++) {
            this->T[i][j] = this->T[i][j] * RO;  //evaporate the pheromones on all possible paths
        }
    }
}


void ACO::calculateVisibilities() {
    std::vector<double> row;
    double visibility, overlapTmp;
    for (int i = 0; i < oligoCount; i++) {
        row.clear();
        for (int j = 0; j < oligoCount; j++) {
            overlapTmp = this->overlap(this->spectrum[i].sequence, this->spectrum[j].sequence, this->oligoLength);
            if (overlapTmp != 0) { //if distance is 0, set visibility as -1
                visibility = 1 / overlapTmp;
            }
            else {
                visibility = -1;
            }

            row.push_back(visibility);
        }
        this->oligoVisibilities.push_back(row);
    }
}


void ACO::calculateOverlaps() {
    std::vector<std::vector<int>> overlaps;
    std::vector<int> row;
    int overlapTmp;
    for (int i = 0; i < oligoCount; i++) {
        row.clear();
        for (int j = 0; j < oligoCount; j++) {
            overlapTmp = this->overlap(this->spectrum[i].sequence, this->spectrum[j].sequence, this->oligoLength);
            row.push_back(overlapTmp);
        }
        this->oligoOverlaps.push_back(row);
    }
}


void ACO::initializePheromones() {
    std::vector<double> row;
    for (int i = 0; i < this->oligoCount; i++) {
        row.clear();
        for (int j = 0; j < this->oligoCount; j++) {
            row.push_back(1.0);
        }
        this->T.push_back(row);
    }
}

void ACO::spawnAnts() {
    for (int i = 0; i < this->oligoCount; i++) {
        Ant ant(this->oligoCount, this->startingOligoIndex);
        this->ants.push_back(ant);
    }
}

std::string ACO::joinSequence(std::vector<Oligo> spectrum) {
    std::string finalSequence = spectrum[0].sequence;
    for (int i=0; i < spectrum.size() - 1; i++) {
        int overlapScore = this->overlap(spectrum[i].sequence, spectrum[i+1].sequence, this->oligoLength);
        if (finalSequence.size() < this->startSeqenceLength) {
            finalSequence.append(spectrum[i+1].sequence.substr(this->oligoLength - overlapScore, overlapScore));
        }
    }
    return finalSequence;
}

int ACO::findStartingOligo(std::string oligo) {
    for (int i=0; i<this->spectrum.size(); i++) {
        if (this->spectrum[i].sequence == oligo) {
            return spectrum[i].index;
        }
    }
    return 0;
}

void ACO::moveAnts() {
    for (int oligo = 0; oligo < this->oligoCount ; oligo++) { //iterate over cities
        for (int l = 0; l < this->antNumber; l++) { //iterate over ants
            this->ants[l].move(T, this->oligoVisibilities); //move each ant
        }
    }
}


void ACO::resetAnts() {
    for (int l = 0; l < this->antNumber; l++) { //iterate over cities
        this->ants[l].reset(this->oligoCount, this->startingOligoIndex);
    }
}

std::string ACO::optimize() {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    int it = 0, resetCounter = 0;

    time_t endwait;
    time_t start = time(NULL);
    time_t sec = 1000;

    endwait = start + sec;

    // Main Loop
    while (start < endwait && resetCounter < this->resetThreshold) {
        auto main_loop_t1 = high_resolution_clock::now();
        std::cout << it << std::endl;
        auto move_loop_t1 = high_resolution_clock::now();
        this->moveAnts();
        auto move_loop_t2 = high_resolution_clock::now();
        duration<double, std::milli> move_loop_time = move_loop_t2 - move_loop_t1;
        std::cout << "move: " << move_loop_time.count() << "ms\n";
        this->updateTrailLevels();
        this->resetAnts();
        this->evaporateTrailLevels();
        // check if the algorithm is providing the same result
        if (this->lastShortestPath == this->shortestPath) {
            resetCounter++;
            std::cout << "Result repeated" << std::endl;
        }
        else { 
            resetCounter = 0;
        }

        it++;
        this->lastShortestPath = this->shortestPath;
        start = time(NULL);
        auto main_loop_t2 = high_resolution_clock::now();
        duration<double, std::milli> main_loop_time = main_loop_t2 - main_loop_t1;
        std::cout << "main: " << main_loop_time.count() << "ms\n";
        std::cout << move_loop_time.count() / main_loop_time.count() * 100 << std::endl;
    }

    std::cout << this->shortestPath.size() << std::endl;

    for (int i=0; i<this->oligoCount; i++){
        std::cout << this->spectrum[this->shortestPath[i]].sequence << " ";
        this->finalSpectrum.push_back(this->spectrum[this->shortestPath[i]]);
    }

    std::string finalSequence = this->joinSequence(this->finalSpectrum);

    return finalSequence;
}