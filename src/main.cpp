#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <map>
#include <chrono>
// #include <windows.h>
#include <set>
#include <fstream>

using namespace std;

double PCFreq = 0.0;
// __int64 CounterStart = 0;

// void StartCounter()
// {
//     LARGE_INTEGER li;
//     if (!QueryPerformanceFrequency(&li))
//         cout << "QueryPerformanceFrequency failed!\n";

//     PCFreq = double(li.QuadPart) / 1000.0;

//     QueryPerformanceCounter(&li);
//     CounterStart = li.QuadPart;
// }
// double GetCounter()
// {
//     LARGE_INTEGER li;
//     QueryPerformanceCounter(&li);
//     return double(li.QuadPart - CounterStart) / PCFreq;
// }

class Constants {
public:
    const float Q = 10;
    const float Ro = 0.8;
    const float alpha = 1;
    const float beta = 1;

    Constants() {};
};

class Oligo {
public:
    int index;
    string sequence;

    Oligo(int index, string sequence) {
        this->index = index;
        this->sequence = sequence;
    }
};

class Ant {
public:
    int current_oligo_index;
    vector <int> path;
    set<int> not_visited;
    double Lk;
    Constants c;

    Ant(/*int current_oligo*/ int n) {
        this->current_oligo_index = 0; //setting the starting point as 1 (for now) ->[spawn_ants]
        this->path.push_back(0);

        //filling the not_visited set, starting from the point after 1 (for now) ->[spawn_ants]
        for (int i = 1; i < n; i++) {
            not_visited.insert(i);
        }
    }

    void reset(int n) {
        this->current_oligo_index = 0;
        this->path.clear();
        this->path.push_back(current_oligo_index);
        this->not_visited.clear();
        for (int i = 1; i < n; i++) {
            not_visited.insert(i);
        }

    }

    void move(vector<vector<double>> T, vector<vector<double>> oligo_visibilities) {
        double random = ((double)rand() / (RAND_MAX));
        double probability;
        double p_n; //numerator probability
        double p_d = 0.0; //denominator probability
        double sum = 0;
        // cout << "random " << random << endl;
        // cout << "not visited size " << not_visited.size() << endl;
        for (set<int>::iterator i = not_visited.begin(); i != not_visited.end(); i++) { //iterate over not visited cities
            p_d = 0.0;
            /*PART 1--------calculate probability for a oligo-----------*/
            p_n = (pow(T[current_oligo_index][*i], c.alpha)) * (pow(oligo_visibilities[current_oligo_index][*i], c.beta)); //numerator probability
            for (set<int>::iterator j = not_visited.begin(); j != not_visited.end(); j++) {
                p_d += (pow(T[current_oligo_index][*j], c.alpha)) * (pow(oligo_visibilities[current_oligo_index][*j], c.beta)); //denominator probability
            }
            probability = p_n / p_d;
            // cout << " probability " << probability << endl;
            // cout << "random " << random << " prob " << probability << endl;

            /*PART 2--------try to select the oligo for which the probability was calculated-----------*/
            random -= probability;
            if (random <= 0) {
                // cout << "random after subtraction " << random << endl;
                current_oligo_index = *i;
                path.push_back(current_oligo_index); //add selected oligo to ant's path
                not_visited.erase(i); //remove selected oligo from not visited
                break;
            }
        }
    }

};


// vector<Oligo> read_spectrum(ifstream infile, string filename) {
//     vector<Oligo> spectrum;
//     string oligo;
//     cout << "Aaa";
//     while (infile >> oligo){
//         spectrum.push_back(Oligo(1, "a"));
//     }

//     return spectrum;
// }


int overlap(std::string a, std::string b, int len) {
    int overlap_score = 0;
    for (int i=1; i<len+1; i++) {
        std::string a_start =  a.substr(0, i);
        std::string b_end = b.substr(len - i, i);

        std::string a_end =  a.substr(len - i, i);
        std::string b_start = b.substr(0, i);

        if (a_start == b_end || a_end == b_start) {
            overlap_score += len - i;  // lower overlap_score == bigger overlap
        }
    }
    
    return overlap_score;
}

double calculate_lk(Ant ant, vector<vector<double>> overlaps, int n) {
    float Lk = 0;
    for (int i = 1; i < n + 1; i++) {
        Lk += overlaps[ant.path[i - 1]][ant.path[i]];
        // cout << point_distances[i-1][i] << " ";
    }
    // cout << endl;
    return Lk;
}

vector<vector<double>> update_trail_levels(Ant ant, double Lk, vector<vector<double>> T, Constants c, int n) {
    vector<vector<double>> pheromone_levels = T;
    for (int i = 1; i < n + 1; i++) { //add the pheromones along the ant's path
        pheromone_levels[ant.path[i - 1]][ant.path[i]] += c.Q / Lk;
    }

    return pheromone_levels;
}


vector<vector<double>> evaporate_trail_levels(vector<vector<double>> T, Constants c, int n) {
    vector<vector<double>> pheromone_levels = T;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            pheromone_levels[i][j] = pheromone_levels[i][j] * c.Ro;  //evaporate the pheromones on all possible paths
        }
    }

    return pheromone_levels;
}


vector<vector<double>> calculate_visibilities(vector<Oligo> oligos, int n, int k) {
    vector<vector<double>> visibilities;
    vector<double> row;
    double visibility, overlap_tmp;
    for (int i = 0; i < n; i++) {
        row.clear();
        for (int j = 0; j < n; j++) {
            overlap_tmp = overlap(oligos[i].sequence, oligos[j].sequence, k);
            if (overlap_tmp != 0) { //if distance is 0, set visibility as -1
                visibility = 1 / overlap_tmp;
            }
            else {
                visibility = -1;
            }

            row.push_back(visibility);
        }
        visibilities.push_back(row);
    }
    return visibilities;
}


vector<vector<double>> calculate_overlaps(vector<Oligo> oligos, int n, int k) {
    vector<vector<double>> overlaps;
    vector<double> row;
    double overlap_tmp;
    for (int i = 0; i < n; i++) {
        row.clear();
        for (int j = 0; j < n; j++) {
            overlap_tmp = overlap(oligos[i].sequence, oligos[j].sequence, k);
            row.push_back(overlap_tmp);
        }
        overlaps.push_back(row);
    }
    return overlaps;
}


vector<vector<double>> initialize_pheromones(int n) {
    vector<vector<double>> pheromones;
    vector<double> row;
    for (int i = 0; i < n; i++) {
        row.clear();
        for (int j = 0; j < n; j++) {
            row.push_back(1.0);
        }
        pheromones.push_back(row);
    }
    return pheromones;
}

vector<Ant> spawn_ants(int n) {
    vector<Ant> ants;
    for (int i = 0; i < n; i++) {
        Ant ant(n);
        ants.push_back(ant);
    }
    return ants;
}


int main(int argc, char * argv[]) {

    srand(time(NULL));
    double minimum = 1000000;
    int iteration_number = 2000;
    vector<int> path;
    string filename = argv[1];
    double Lk;
    Constants c;
    vector<int> shortest_path;

    // read parameters
    ifstream infile(filename);
    int n, k, index = 0;
    int ant_number = n;
    infile >> n;
    infile >> k;
    
    // read oligos
    vector<Oligo> spectrum;
    string sequence;
    while (infile >> sequence){
        spectrum.push_back(Oligo(index, sequence));
        index++;
    }

    vector<vector<double>> oligo_overlaps = calculate_overlaps(spectrum, n, k); //initialize oligo overlap table
    vector<vector<double>> oligo_visibilities = calculate_visibilities(spectrum, n, k); //initialize oligo visibilities table
    vector<vector<double>> T = initialize_pheromones(n); //initialize pheromones table

    /* ------test oligo visibilities------ */
    // for (int i=0; i<n; i++){
    //     for (int j=0; j<n; j++){
    //         cout << oligo_visibilities[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    /* ------test oligo distances------ */
    // for (int i=0; i<n; i++){
    //     for (int j=0; j<n; j++){
    //         cout << oligo_distances[i][j] << " ";
    //     }
    //     cout << endl;
    // }


    /* ------test pheromones------ */
    // for (int i=0; i<n; i++){
    //     for (int j=0; j<n; j++){
    //         cout << T[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    vector<Ant> ants = spawn_ants(n); //spawn all ants

    //MAIN LOOP
    for (int t = 0; t < iteration_number; t++) {
        cout << "Iteration: " << t << endl;

        // MOVE LOOP
        for (int oligo = 0; oligo < n; oligo++) { //iterate over cities
            for (int l = 0; l < ant_number; l++) { //iterate over ants
                ants[l].move(T, oligo_visibilities); //move each ant
            }
        }

        // RESET AND UPDATE LOOP
        for (int l = 0; l < ant_number; l++) { //iterate over ants
            ants[l].path.push_back(0); //return to spawn oligo
            Lk = calculate_lk(ants[l], oligo_overlaps, n); //calculate the total distance for ant
            T = update_trail_levels(ants[l], Lk, T, c, n); //add the pheromones along the ant's path
            if (Lk < minimum) {
                minimum = Lk; //check if the path is the shortest yet
                // shortest_path = ants[l].path;
            }
            ants[l].reset(n);

        }
        T = evaporate_trail_levels(T, c, n); //evaporate all pheromones
        // cout << "Shortest path: ";
        // for (int i=0; i<n; i++){
        //     cout << shortest_path[i] << " ";
        // }
        // cout << endl;
        /* ------test pheromones------ */
        // for (int i = 0; i < n; i++) {
        //     for (int j = 0; j < n; j++) {
        //         cout << T[i][j] << " ";
        //     }
        //     cout << endl;
        // }
    }
    cout << "Minimum: " << minimum << endl;
}


