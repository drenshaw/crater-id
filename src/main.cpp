#include <iostream>
#include <Eigen/Dense>
#include "structs.h"
#include "crater-id.h"
#include "combinatorics.h"
#include "io.h"
#include "conics.h"
// #include "gnuplot-iostream.h"
 
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
    Conic conicD(12, 8, -500, -20, 0);
    // std::tuple<Eigen::Vector3d, Eigen::Vector3d> gh_ij, gh_jk, gh_ki;
    // Eigen::Vector3d lij, ljk, lki;
    // double invA, invB, invC;
    // // conicA.ConicIntersectionLines(conicB, gh);
    // if(!IntersectionLines(conicA.getLocus(), conicB.getLocus(), gh_ij)) {std::cout<<"1\n";return 1;}
    // if(!IntersectionLines(conicB.getLocus(), conicC.getLocus(), gh_jk)) {std::cout<<"2\n";return 2;}
    // if(!IntersectionLines(conicC.getLocus(), conicA.getLocus(), gh_ki)) {std::cout<<"3\n";return 3;}
    // if(!ChooseIntersection(gh_ij, conicA.getCenter(), conicB.getCenter(), lij)) {std::cout<<"4\n";return 4;}
    // if(!ChooseIntersection(gh_jk, conicB.getCenter(), conicC.getCenter(), ljk)) {std::cout<<"5\n";return 5;}
    // if(!ChooseIntersection(gh_ki, conicC.getCenter(), conicA.getCenter(), lki)) {std::cout<<"6\n";return 6;}
    // if(!computeInvariant(lij, lki, conicA.getLocus(), invA)) {std::cout<<"7\n";return 7;}
    // if(!computeInvariant(ljk, lij, conicB.getLocus(), invB)) {std::cout<<"8\n";return 8;}
    // if(!computeInvariant(lki, ljk, conicC.getLocus(), invC)) {std::cout<<"9\n";return 9;}
    
    // std::cout << "Invariants: " << invA << ", " << invB << ", " << invC << std::endl;

    std::vector<double> invariantsABC, invariantsBCD, invariantsCDA, invariantsDAB;
    if(!computeCraterTriadInvariants(conicA, conicB, conicC, invariantsABC)) {
        std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
    }
    if(!computeCraterTriadInvariants(conicB, conicC, conicD, invariantsBCD)) {
        std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
    }
    if(!computeCraterTriadInvariants(conicC, conicD, conicA, invariantsCDA)) {
        std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
    }
    if(!computeCraterTriadInvariants(conicD, conicA, conicB, invariantsDAB)) {
        std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
    }
    printVector(invariantsABC, "Invariants: ");
    printVector(invariantsBCD, "Invariants: ");
    printVector(invariantsCDA, "Invariants: ");
    printVector(invariantsDAB, "Invariants: ");
    return 0;
}

