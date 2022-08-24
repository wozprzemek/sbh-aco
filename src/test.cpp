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


int main() {
    double a = 0;
    for (int i=0; i < 500; i++) {
        for (int j=0; j< 500; j++) {
            for (int k=0; k < 500; k++) {
                if (k % 2 == 0) {
                    a += 1.001;
                }
                else {
                    a--;
                }
            }
        }
    }
    std::cout << a;
}