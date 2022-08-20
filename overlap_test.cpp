#include <iostream>
#include <algorithm> 

int overlap(std::string a, std::string b, int len) {
    int overlap_val = 0;
    for (int i=1; i<len+1; i++) {
        std::string a_start =  a.substr(0, i);
        std::string b_end = b.substr(len - i, i);

        std::string a_end =  a.substr(len - i, i);
        std::string b_start = b.substr(0, i);

        if (a_start == b_end || a_end == b_start) {
            overlap_val += len - i;
        }
    }
    
    return overlap_val;
}

int main() {

    std::string oligo1 = "GTAC";
    std::string oligo2 = "CATG";
    
    int overlap_val = overlap(oligo1, oligo2, 4);
    std::cout << overlap_val << std::endl;
}