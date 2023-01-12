#include <iostream>
#include "structs.h"
#include "crater-id.h"
#include "combinatorics.h"
#include "io.h"
#include "conics.h"
#include <Eigen/Dense>
 
int main() {
    // std::string fname;
    // const std::string degrees = "°";
    // // cout<<"Enter crater file name: ";
    // // cin>>fname;
    // // fname = "/home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv";
    // fname = "/home/dqr0509/data/craters/lunar_craters.csv";

    // std::vector<lunar_crater> craters;
    // runCraterReader(fname, craters);
    
    // std::vector<std::tuple<uint, uint>> valids;
    // specialCombination(craters, valids, 10.);
    
    // std::vector<std::tuple<uint, uint, uint>> triads;
    // formTriads(valids, triads);
    // std::cout << triads.size() << " valid triads found." << std::endl;
    // // uint idx = 0;
    // // for(const auto& [i, j, k]: triads) {
    // //     std::cout << "IDX: " << idx++ << " | " << i << ", " << j << ", " << k << std::endl;
    // // }

    // uint count = 0;
    // for(const auto& [i, j, k] : triads) {
    //     std::cout << "-" << count++ << "->"
    //               << " "  << i
    //               << ", " << j
    //               << ", " << k
    //               << std::endl;
    // }
    
    Conic conicA(10, 7, -100, -50, 0);
    Conic conicB(15, 12, 100, 200, 0);
    Conic conicC(12, 8, 50, -200, 0);
    // std::cout << "Conic A: " << conicA.getLocus() << std::endl;
    // std::cout << "Conic B: " << conicB.getLocus() << std::endl;
    std::tuple<Eigen::Vector3f, Eigen::Vector3f> gh_ij, gh_jk, gh_ki;
    Eigen::Vector3f lij, ljk, lki;
    float invA, invB, invC;
    // conicA.ConicIntersectionLines(conicB, gh);
    if(!IntersectionLines(conicA.getLocus(), conicB.getLocus(), gh_ij)) {std::cout<<"1\n";return 1;}
    std::cout << std::endl << "*****************" << std::endl;
    if(!IntersectionLines(conicB.getLocus(), conicC.getLocus(), gh_jk)) {std::cout<<"2\n";return 2;}
    std::cout << std::endl << "*****************" << std::endl;
    if(!IntersectionLines(conicC.getLocus(), conicA.getLocus(), gh_ki)) {std::cout<<"3\n";return 3;}
    std::cout << std::endl << "*****************" << std::endl;
    if(!ChooseIntersection(gh_ij, conicA.getCenter(), conicB.getCenter(), lij)) {std::cout<<"4\n";return 4;}
    std::cout << std::endl << "++++++++++++++" << std::endl;
    if(!ChooseIntersection(gh_jk, conicB.getCenter(), conicC.getCenter(), ljk)) {std::cout<<"5\n";return 5;}
    std::cout << std::endl << "++++++++++++++" << std::endl;
    if(!ChooseIntersection(gh_ki, conicC.getCenter(), conicA.getCenter(), lki)) {std::cout<<"6\n";return 6;}
    std::cout << std::endl << "++++++++++++++" << std::endl;
    std::cout << "Line_ij: \n" << lij << std::endl;
    std::cout << "Line_jk: \n" << ljk << std::endl;
    std::cout << "Line_ki: \n" << lki << std::endl;
    // if(!computeInvariant(lij, lki, conicA.getLocus(), invA)) {
    //     return 7;
    // }
    // if(!computeInvariant(ljk, lij, conicB.getLocus(), invB)) {
    //     return 8;
    // }
    if(!computeInvariant(lki, ljk, conicC.getLocus(), invC)) {
        return 9;
    }
    
    return 0;
}

