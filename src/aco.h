#include <vector>
#include <set>
#include <cstdlib>
#include "constants.h"
#include <string>
#include <math.h>

class Ant {
public:
    int currentOligoIndex;
    std::vector <int> path;
    std::set<int> notVisited;
    double Lk, alpha, beta, rho;

    Ant(int oligoCount, int startingIndex, double alpha, double beta, double rho);

    void reset(int oligoCount, int startingIndex);
    void move(std::vector<std::vector<double>> T, std::vector<std::vector<double>> oligoVisibilities);

};

class Oligo {
public:
    int index;
    std::string sequence;

    Oligo(int index, std::string sequence);
};


class ACO {
public:
    int startSeqenceLength, oligoLength, oligoCount, antNumber, startingOligoIndex, resetThreshold=20, minimum;
    double Lk, alpha, beta, rho;
    std::string startSequence, startingOligo;
    std::vector<Oligo> spectrum, finalSpectrum;
    std::vector<int> shortestPath, lastShortestPath;
    std::vector<std::vector<int>> oligoOverlaps;
    std::vector<std::vector<double>> oligoVisibilities;
    std::vector<std::vector<double>> T;
    std::vector<Ant> ants;
    

	ACO(int startSeqenceLength, int oligoLength, std::string startSequence, std::vector<Oligo> spectrum, int antNumber, double alpha, double beta, double rho);

    void init();
    std::string optimize();
    int overlap(std::string a, std::string b, int len);
    double calculateLk(Ant ant);

    void calculateVisibilities();
    void calculateOverlaps();
    void initializePheromones();
    void spawnAnts();
    void moveAnts();
    void resetAnts();
    void updateTrailLevels();
    void evaporateTrailLevels();

    std::string joinSequence(std::vector<Oligo> finalSpectrum);
    int findStartingOligo(std::string oligo);
};