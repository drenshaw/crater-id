#include <iostream>
#include "structs.h"
#include "crater-id.h"
#include "combinatorics.h"
#include "io.h"
#include "conics.h"
 
int main() {
    std::string fname;
    const std::string degrees = "°";
    // cout<<"Enter crater file name: ";
    // cin>>fname;
    // fname = "/home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv";
    // fname = "/home/dqr0509/data/craters/lunar_craters.csv";

    // std::vector<lunar_crater> craters;
    // runCraterReader(fname, craters);
    
    // std::vector<std::tuple<uint, uint>> valids;
    // specialCombination(craters, valids, 15.);
    
    // std::vector<std::tuple<uint, uint, uint>> triads;
    // formTriads(valids, triads);
    // std::cout << triads.size() << " valid triads found." << std::endl;
    
    Conic conic(50, 40, 10, 10, 0);
    std::vector<float> coeff;
    std::vector<float> geom = conic.getGeom();
    printVector(geom, "Geom: ");
    coeff = conic.Geom2Implicit();
    printVector(coeff, "Coeff: ");
    // uint idx = 0;
    // for(const auto& [i, j, k]: triads) {
    //     std::cout << "IDX: " << idx++ << " | " << i << ", " << j << ", " << k << std::endl;
    // }

    // uint count = 0;
    // for(const auto& [i, j, k] : triads) {
    //     std::cout << "-" << count++ << "->"
    //               << " " << craters[i]
    //               << ", " << craters[j]
    //               << ", " << craters[k]
    //               << std::endl;
    // }
    return 0;
}
