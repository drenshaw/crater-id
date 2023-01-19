#pragma once

#include <iostream>
#include "structs.h"
#include "io.h"

std::ostream& operator<<(std::ostream& os, const lunar_crater& crater);
std::ostream& operator<<(std::ostream& os, const Point& point);

// void printLunarCratersInfo(const lunar_crater);
// bool readLunarCraterEntry(const std::string, const char);
// template <typename T>
// void Permutation(std::vector<T>);
// template <typename T>
// std::vector<std::vector<T>> Combination(const std::vector<T>, const uint);
// template <typename T, typename N>
// std::vector<std::tuple<T, T>> specialCombination(const std::vector<T>, const N =30.);
// template <typename T, typename N>
// void specialCombination(const std::vector<T>, 
//                         const N,
//                         std::vector<std::tuple<uint, uint>>&);
// std::vector<std::tuple<uint, uint, uint>> formTriad(const std::vector<std::tuple<uint, uint>>);
// std::ostream& operator<<(std::ostream&, const lunar_crater&);
// std::ostream& operator<<(std::ostream&, const Point&);
// template <typename T>
// T deg2rad(const T);
// template <typename T>
// T rad2deg(const T);
// template <typename T>
// Point latlon2unitVector(const T lat, const T lon);
// template <typename T>
// Point latlon2unitVector(const T crater);
// template <typename R, typename T>
// Point LLHtoECEF(const R, const T);
// double vectorNorm(const Point pt);
// void normalizeVector(Point&);
// double dot(const Point, const Point);
// double angularDistance(const Point, const Point);
// double angularPseudodistance(const Point, const Point);
// template <typename T>
// double latlon_dist(const T, const T, const T, const T);
// template <typename T>
// double latlon_dist(const T, const T);