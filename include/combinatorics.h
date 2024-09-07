#ifndef COMBINATORICS_H
#define COMBINATORICS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <iterator>

#include "structs.h"
#include "crater-id.h"
#include "vector_math.h"

// template <typename T, typename N>
// std::vector<std::tuple<T, T>> specialCombination(const std::vector<T>, 
//                                                  const N);
template <typename T, typename N>
void specialCombination(const std::vector<T>,
                        std::vector<std::tuple<uint, uint>>&, 
                        const N=10.);
template <typename T>
std::vector<std::vector<T>> Combination(const std::vector<T>, const uint);
template <typename T>
void Permutation(std::vector<T> v);

void formTriads(const std::vector<std::tuple<uint, uint>>,
                std::vector<std::tuple<uint, uint, uint>>&);

/**** Template definitions ****/
template <typename T, typename N>
void specialCombination(const std::vector<T> choices, 
                        std::vector<std::tuple<uint, uint>>& valid_craters,
                        const N max_angle_deg) {
    const uint n = choices.size();
    // const double max_angle_deg = 30.;
    const N max_angle_rad = deg2rad(max_angle_deg);
    // trig function are expensive in loops, so use the max dot product
    // angle = acos(dot(a,b)) == cos(angle) = dot(a,b)
    const N min_dot_prod = cos(max_angle_rad);

    Eigen::Vector3d pt1, pt2;
    // std::vector<std::tuple<T, T>> valid_craters;
    std::string lat1, lat2, lon1, lon2;
    T current_choice, next_choice;

    // uint count = 0;
    for(size_t i=0; i<n-1; i++) {
        current_choice = choices[i];
        pt1 = latlon2bearing(current_choice);
        for(size_t j=i+1; j<n; j++) {
            next_choice = choices[j];
            pt2 = latlon2bearing(next_choice);
            // if(angularDistance(pt1, pt2) < max_angle_rad) {
            if(angularPseudoDistance(pt1, pt2) > min_dot_prod) {
                // count++;
                valid_craters.push_back({i, j});
            }
        }
    }
    std::cout << valid_craters.size() << " valid pairs." << std::endl;
}

template <typename T, typename N>
std::vector<std::tuple<T, T>> specialCombination(const std::vector<T> choices, 
                                                 const N max_angle_deg) {
    const uint n = choices.size();
    // const double max_angle_deg = 30.;
    const N max_angle_rad = deg2rad(max_angle_deg);
    // trig function are expensive in loops, so use the max dot product
    // angle = acos(dot(a,b)) == cos(angle) = dot(a,b)
    const N min_dot_prod = cos(max_angle_rad);

    Eigen::Vector3d pt1, pt2;
    std::vector<std::tuple<T, T>> valid_craters;
    std::string lat1, lat2, lon1, lon2;
    T current_choice, next_choice;

    // uint count = 0;
    for(size_t i=0; i<n-1; i++) {
        current_choice = choices[i];
        pt1 = latlon2bearing(current_choice);
        for(size_t j=i+1; j<n; j++) {
            next_choice = choices[j];
            pt2 = latlon2bearing(next_choice);
            // if(angularDistance(pt1, pt2) < max_angle_rad) {
            if(angularPseudoDistance(pt1, pt2) > min_dot_prod) {
                // count++;
                valid_craters.push_back({current_choice, next_choice});
            }
        }
    }
    std::cout << "Valid combinations: " << valid_craters.size() << std::endl;
    return valid_craters;
}

void print_triads(const std::vector<std::tuple<uint, uint, uint>> triads, 
                  const std::vector<lunar_crater> craters,
                  const uint max_iter=10);

#endif