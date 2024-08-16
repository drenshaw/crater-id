#ifndef VECTOR_MATH
#define VECTOR_MATH

#include <math.h>
#include <vector>
#include <array>
// #include <experimental/array>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include "crater-id.h"

#define EPS (10 * std::numeric_limits<T>::epsilon())

double getCofactor(const Eigen::MatrixXd& matrix, int p, int q);
Eigen::MatrixXd getCofactorMatrix(const Eigen::MatrixXd& matrix);
Eigen::MatrixXd getMatrixAdjugate(const Eigen::MatrixXd&);
Eigen::Matrix3d get3x3SymmetricMatrixAdjugate(const Eigen::Matrix3d&);

// TODO: be careful using templates here: if an "int" is passed, we get an "int" back
template <typename T>
double deg2rad(const T deg);
template <typename T>
double rad2deg(const T rad);
template <typename T>
Eigen::Vector3d latlon2bearing(const T lat, const T lon);
template <typename T>
Eigen::Vector3d latlon2bearing(const T crater);
// TODO: Why did I have these overloaded operators?
// template <typename T>
// void operator /(Eigen::Vector3d& vec, T divisor);
// template <typename T>
// void operator *(Eigen::Vector3d& vec, T scalar);
template <typename T>
Eigen::Vector3d llarToWorld(const T crater, double alt, double radius=1.0);
template <typename R, typename T>
Eigen::Vector3d LLHtoECEF(const R crater, const T alt);
// double vectorNorm(const Eigen::Vector3d& pt);
// void normalizeVector(Eigen::Vector3d& pt);
// template<typename Iter_T>
// long double vectorNorm(Iter_T first, Iter_T last);
// template<typename T>
// T vectorNorm(std::vector<T> vec);
double vdot(const Eigen::Vector3d& pt1, const Eigen::Vector3d& pt2);
double angularPseudoDistance(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2);
double angularDistance(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2);
template <typename T>
double latlon_dist(const T crater1, const T crater2);
template <typename T>
double latlon_dist(const T lat1, const T lon1, const T lat2, const T lon2);
template <typename T>
void makeUnique(T& vec);
template <typename T>
std::vector<uint> getRange(std::vector<T> vec);
bool normalizeDeterminant(Eigen::MatrixXd& mtx);
Eigen::Matrix3d crossMatrix(const Eigen::Vector3d&);
void normalizeVector(Eigen::Vector3d&);
void normalizeVector(const Eigen::Vector3d&, Eigen::Vector3d&);
Eigen::Vector3d normalizeVector(const Eigen::Vector3d&);
template <typename T>
T vectorNorm(const Eigen::Vector3d&);
void GetNorthPoleUnitVector(Eigen::Vector3d&);
Eigen::Vector3d GetNorthPoleUnitVector();

/**** Template definitions ****/

template <typename T, typename A>
int arg_max(std::vector<T, A> const& vec) {
  return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
}

template <typename T, typename A>
int arg_min(std::vector<T, A> const& vec) {
  return static_cast<int>(std::distance(vec.begin(), min_element(vec.begin(), vec.end())));
}

template <typename T, size_t SIZE>
int arg_max(std::array<T, SIZE> const& arr) {
  return static_cast<int>(std::distance(arr.begin(), std::max_element(arr.begin(), arr.end())));
}

template <typename T, size_t SIZE>
int arg_min(std::array<T, SIZE> const& arr) {
  return static_cast<int>(std::distance(arr.begin(), std::min_element(arr.begin(), arr.end())));
}

template <typename T>
double deg2rad(const T deg) {
  return static_cast<double>(deg) * M_PI/180.;
}
template <typename T>
double rad2deg(const T rad) {
  return static_cast<float>(rad) * 180. / M_PI;
}

template <typename T>
Eigen::Vector3d latlon2bearing(const T lat, const T lon) {
    double lat_rad = deg2rad(lat);
    double lon_rad = deg2rad(lon);

    Eigen::Vector3d point;
    point << cos(lat_rad) * cos(lon_rad),
             cos(lat_rad) * sin(lon_rad),
             sin(lat_rad);

    return point;
}

// TODO: we don't use this currently
template <typename T>
Eigen::Vector3d latlon2bearing(const T crater) {
    return latlon2bearing(crater.lat, crater.lon);
}

template <typename Derived, int size>
Derived vectorNorm(const Eigen::Vector<Derived, size>& vec) {
  return vec/vec.norm();
}

template <typename Iter_T>
long double vectorNorm(Iter_T first, Iter_T last) {
  return sqrt(std::inner_product(first, last, first, 0.0L));
}

template <typename T>
T vectorNorm(const std::vector<T> vec) {
    return vectorNorm(vec.begin(), vec.end());
}

template <typename T, size_t SIZE>
T vectorNorm(const std::array<T, SIZE> vec) {
    return vectorNorm(vec.begin(), vec.end());
}

template <typename T, size_t SIZE>
void copy_vec2array(const std::vector<T> vec, std::array<T, SIZE>& arr) {
  std::copy_n(std::make_move_iterator(vec.begin()), SIZE, arr.begin());
}

#endif