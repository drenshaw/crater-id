#pragma once

#include <math.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include "crater-id.h"

template <typename T>
T deg2rad(const T deg);
template <typename T>
T rad2deg(const T rad);
template <typename T>
Eigen::Vector3f latlon2unitVector(const T lat, const T lon);
template <typename T>
Eigen::Vector3f latlon2unitVector(const T crater);
template <typename T>
Eigen::Vector3f llarToWorld(const T crater, float alt, float rad=1.0);
template <typename R, typename T>
Eigen::Vector3f LLHtoECEF(const R crater, const T alt);
// float vectorNorm(const Eigen::Vector3f& pt);
// void normalizeVector(Eigen::Vector3f& pt);
// template<typename Iter_T>
// long float vectorNorm(Iter_T first, Iter_T last);
// template<typename T>
// T vectorNorm(std::vector<T> vec);
float vdot(const Eigen::Vector3f& pt1, const Eigen::Vector3f& pt2);
float angularPseudodistance(const Eigen::Vector3f& point1, const Eigen::Vector3f& point2);
float angularDistance(const Eigen::Vector3f& point1, const Eigen::Vector3f& point2);
template <typename T>
float latlon_dist(const T crater1, const T crater2);
template <typename T>
float latlon_dist(const T lat1, const T lon1, const T lat2, const T lon2);
template <typename T>
void makeUnique(T& vec);
template <typename T>
std::vector<uint> getRange(std::vector<T> vec);
Eigen::Matrix3f normalizeDeterminant(const Eigen::Matrix3f& mtx);
Eigen::Matrix3f crossMatrix(const Eigen::Vector3f&);

/**** Template definitions ****/

template <typename T, typename A>
int arg_max(std::vector<T, A> const& vec) {
  return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
}

template <typename T, typename A>
int arg_min(std::vector<T, A> const& vec) {
  return static_cast<int>(std::distance(vec.begin(), min_element(vec.begin(), vec.end())));
}

template <typename T>
T deg2rad(const T deg) {
    return deg * M_PI/180.;
}
template <typename T>
T rad2deg(const T rad) {
    return rad * 180. / M_PI;
}

template <typename T>
Eigen::Vector3f latlon2unitVector(const T lat, const T lon) {
    float lat_rad = deg2rad(lat);
    float lon_rad = deg2rad(lon);

    Eigen::Vector3f point;
    point << cos(lat_rad) * cos(lon_rad),
             cos(lat_rad) * sin(lon_rad),
             sin(lat_rad);

    return point;
}

template <typename T>
Eigen::Vector3f latlon2unitVector(const T crater) {
    return latlon2unitVector(crater.lat, crater.lon);
}

template<typename Iter_T>
long double vectorNorm(Iter_T first, Iter_T last) {
  return sqrt(inner_product(first, last, first, 0.0L));
}

template<typename T>
T vectorNorm(const std::vector<T> vec) {
    return vectorNorm(vec.begin(), vec.end());
}
