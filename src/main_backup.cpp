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
using namespace std::chrono_literals;

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

    Ant(int oligo_count, int starting_index) {
        this->current_oligo_index = starting_index; //setting the starting point
        this->path.push_back(starting_index);

        //filling the not_visited set
        for (int i = 0; i < oligo_count; i++) {
            if (i != starting_index) {
                not_visited.insert(i);
            }
        }
    }

    void reset(int oligo_count, int starting_index) {
        this->current_oligo_index = starting_index;
        this->path.clear();
        this->path.push_back(current_oligo_index);
        this->not_visited.clear();
        for (int i = 0; i < oligo_count; i++) {
            if (i != starting_index) {
                not_visited.insert(i);
            }
        }

    }

    void move(vector<vector<double>> T, vector<vector<double>> oligo_visibilities) {
        double random = ((double)rand() / (RAND_MAX));
        double probability;
        double p_n; //numerator probability
        double p_d = 0.0; //denominator probability
        double sum = 0;

        for (set<int>::iterator j = not_visited.begin(); j != not_visited.end(); j++) {
            p_d += (pow(T[current_oligo_index][*j], c.alpha)) * (pow(oligo_visibilities[current_oligo_index][*j], c.beta)); //denominator probability
        }

        for (set<int>::iterator i = not_visited.begin(); i != not_visited.end(); i++) { //iterate over not visited cities
            /*PART 1--------calculate probability for a oligo-----------*/
            p_n = (pow(T[current_oligo_index][*i], c.alpha)) * (pow(oligo_visibilities[current_oligo_index][*i], c.beta)); //numerator probability
            probability = p_n / p_d;
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

int overlap(std::string a, std::string b, int len) {
    int overlap_score = len;
    // cout << "FUN overlap: " << a << ", " << b << " " << len << endl;
    for (int i=1; i<len+1; i++) {
        // std::string a_start =  a.substr(0, i);
        // std::string b_end = b.substr(len - i, i);

        std::string a_end =  a.substr(len - i, i);
        std::string b_start = b.substr(0, i);

        if (a_end == b_start) {
            overlap_score = len - i; // lower overlap_score == bigger overlap
        }
    }
    
    return overlap_score;
}

double calculate_lk(Ant ant, vector<vector<int>> overlaps, int oligo_count) {
    float Lk = 0;
    for (int i = 1; i < oligo_count; i++) {
        Lk += overlaps[ant.path[i - 1]][ant.path[i]];
        // cout << point_distances[i-1][i] << " ";
    }
    // cout << endl;
    return Lk;
}

vector<vector<double>> update_trail_levels(Ant ant, double Lk, vector<vector<double>> T, Constants c, int oligo_count) {
    vector<vector<double>> pheromone_levels = T;
    for (int i = 1; i < oligo_count; i++) { //add the pheromones along the ant's path
        pheromone_levels[ant.path[i - 1]][ant.path[i]] += c.Q / Lk; 
    }

    return pheromone_levels;
}


vector<vector<double>> evaporate_trail_levels(vector<vector<double>> T, Constants c, int oligo_count) {
    vector<vector<double>> pheromone_levels = T;

    for (int i = 0; i < oligo_count; i++) {
        for (int j = 0; j < oligo_count; j++) {
            pheromone_levels[i][j] = pheromone_levels[i][j] * c.Ro;  //evaporate the pheromones on all possible paths
        }
    }

    return pheromone_levels;
}


vector<vector<double>> calculate_visibilities(vector<Oligo> oligos, int oligo_count, int k) {
    vector<vector<double>> visibilities;
    vector<double> row;
    double visibility, overlap_tmp;
    for (int i = 0; i < oligo_count; i++) {
        row.clear();
        for (int j = 0; j < oligo_count; j++) {
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


vector<vector<int>> calculate_overlaps(vector<Oligo> oligos, int oligo_count, int k) {
    // cout << "FUN calculate_overlaps" << endl;
    vector<vector<int>> overlaps;
    vector<int> row;
    int overlap_tmp;
    for (int i = 0; i < oligo_count; i++) {
        row.clear();
        for (int j = 0; j < oligo_count; j++) {
            // cout << i << " " << j << endl;
            overlap_tmp = overlap(oligos[i].sequence, oligos[j].sequence, k);
            // cout << oligos[i].sequence << " " << oligos[j].sequence << " " << overlap_tmp << endl;
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

vector<Ant> spawn_ants(int oligo_count, int starting_oligo_index) {
    vector<Ant> ants;
    for (int i = 0; i < oligo_count; i++) {
        Ant ant(oligo_count, starting_oligo_index);
        ants.push_back(ant);
    }
    return ants;
}

string join_sequence(vector<Oligo> final_spectrum, int len, int n) {
    string final_sequence = final_spectrum[0].sequence;
    for (int i=0; i < final_spectrum.size() - 1; i++) {
        int overlap_score = overlap(final_spectrum[i].sequence, final_spectrum[i+1].sequence, len);
        // final_sequence.append(final_spectrum[i].sequence);
        if (final_sequence.size() < n) {
            final_sequence.append(final_spectrum[i+1].sequence.substr(len - overlap_score, overlap_score));
        }
    }

    return final_sequence;
}

int find_starting_oligo(string oligo,vector<Oligo> spectrum) {
    for (int i=0; i<spectrum.size(); i++) {
        if (spectrum[i].sequence == oligo) {
            return spectrum[i].index;
        }
    }
    return 0;
}

// double calculate_accuracy(string start_sequence, string final_sequence) {
//     for 
// }

int main(int argc, char * argv[]) {

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    srand(time(NULL));
    double minimum = 1000000;
    int iteration_number = 100;
    vector<int> path;
    string filename = argv[1];
    double Lk;
    Constants c;
    vector<int> shortest_path;
    vector<int> last_shortest_path;
    string start_sequence;
    string starting_oligo;
    int starting_oligo_index;
    
    // read parameters
    ifstream infile(filename);
    ofstream outfile("../src/output.txt");

    int n, k, index = 0;
    infile >> n;
    infile >> k;
    infile >> start_sequence;


    // read oligos
    vector<Oligo> spectrum;
    vector<Oligo> final_spectrum;

    string sequence;
    while (infile >> sequence){
        spectrum.push_back(Oligo(index, sequence));
        index++;
    }

    int oligo_count = spectrum.size();
    int ant_number = oligo_count;

    infile.close();

    starting_oligo = start_sequence.substr(0, k);
    starting_oligo_index = find_starting_oligo(starting_oligo, spectrum);

    cout << "INITIALIZATION START" << endl;
    vector<vector<int>> oligo_overlaps = calculate_overlaps(spectrum, oligo_count, k); //initialize oligo overlap table
    vector<vector<double>> oligo_visibilities = calculate_visibilities(spectrum, oligo_count, k); //initialize oligo visibilities table
    vector<vector<double>> T = initialize_pheromones(n); //initialize pheromones table

    cout << "INITIALIZATION SUCCESSFUL" << endl;
    // cout << "n: " << n << endl;
    // cout << "spectrum size: " << spectrum.size() << endl;

    /* ------test oligo visibilities------ */
    // for (int i=0; i<oligo_count; i++){
    //     for (int j=0; j<oligo_count; j++){
    //         cout << oligo_visibilities[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    /* ------test oligo overlaps------ */
    // cout << oligo_overlaps[0].size() << endl;

    // for (int i=0; i<oligo_count; i++){
    //     for (int j=0; j<oligo_count; j++){
    //         cout << oligo_overlaps[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    // /* ------test pheromones------ */
    // for (int i=0; i<oligo_count; i++){
    //     for (int j=0; j<oligo_count; j++){
    //         cout << T[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    vector<Ant> ants = spawn_ants(ant_number, starting_oligo_index); //spawn all ants

    int t = 0;
    int reset_counter = 0;
    int reset_threshold = 10;
    // cout << "ANTS SPAWNED" << endl;

    //MAIN LOOP
    // for (int t = 0; t < iteration_number; t++) {
    while (t < iteration_number && reset_counter < reset_threshold) {
        auto main_loop_t1= high_resolution_clock::now();

        cout << "ITERATION " << t << endl;

        auto move_loop_t1 = high_resolution_clock::now();

        // MOVE LOOP
        for (int oligo = 0; oligo < oligo_count; oligo++) { //iterate over cities
            for (int l = 0; l < ant_number; l++) { //iterate over ants
                ants[l].move(T, oligo_visibilities); //move each ant
            }
        }

        cout << "MOVE LOOP" << endl;
        auto move_loop_t2 = high_resolution_clock::now();
        duration<double, std::milli> move_loop_time = move_loop_t2 - move_loop_t1;
        cout << move_loop_time.count() << "ms\n";

        // cout << "   ANTS MOVED" << endl;

        auto update_loop_t1 = high_resolution_clock::now();

        // RESET AND UPDATE LOOP
        for (int l = 0; l < ant_number; l++) { //iterate over ants
            Lk = calculate_lk(ants[l], oligo_overlaps, oligo_count); //calculate the total distance for ant
            T = update_trail_levels(ants[l], Lk, T, c, oligo_count); //add the pheromones along the ant's path
            // cout << Lk << " ";
            if (Lk < minimum) {
                minimum = Lk; //check if the path is the shortest yet
                shortest_path = ants[l].path;
            }

            ants[l].reset(oligo_count, starting_oligo_index);

        }

        cout << "UPDATE LOOP" << endl;
        auto update_loop_t2 = high_resolution_clock::now();
        duration<double, std::milli> update_loop_time = update_loop_t2 - update_loop_t1;
        cout << update_loop_time.count() << "ms\n";

        // cout << "   TRAILS UPDATED" << endl;

        auto evaporate_t1 = high_resolution_clock::now();

        T = evaporate_trail_levels(T, c, oligo_count); //evaporate all pheromones
        if (last_shortest_path == shortest_path) {
            cout << "Same!" << endl;
            reset_counter++;
        }
        else { 
            reset_counter = 0;
        }

        cout << "EVAPORATE" << endl;
        auto evaporate_t2 = high_resolution_clock::now();
        duration<double, std::milli> evaporate_time = evaporate_t2 - evaporate_t1;
        cout << evaporate_time.count() << "ms\n";

        t++;
        last_shortest_path = shortest_path;
        
        cout << "MAIN LOOP" << endl;
        auto main_loop_t2 = high_resolution_clock::now();
        duration<double, std::milli> main_loop_time = main_loop_t2 - main_loop_t1;
        cout << main_loop_time.count() << "ms\n\n";
    }
    // cout << "Minimum: " << minimum << endl;
    cout << "--Final spectrum:  ";
    for (int i=0; i<oligo_count; i++){
        cout << spectrum[shortest_path[i]].sequence << " ";
        final_spectrum.push_back(spectrum[shortest_path[i]]);
    }

    string final_sequence = join_sequence(final_spectrum, k, n);

    cout << endl << endl;
    cout << start_sequence << endl;
    cout << final_sequence << endl << endl;
    outfile << start_sequence << endl << final_sequence;
    outfile.close();
}


