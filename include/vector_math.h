#ifndef VECTOR_MATH
#define VECTOR_MATH

#include <math.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include "crater-id.h"

#define EPS (10 * std::numeric_limits<T>::epsilon())

template <typename Derived, int size>
Derived getCofactor(const Eigen::Matrix<Derived, size, size>& matrix, size_t cf_row, size_t cf_col);
template <typename Derived, int size>
Eigen::Matrix<Derived, size, size> getCofactorMatrix(const Eigen::Matrix<Derived, size, size>& matrix);
template <typename Derived, int size>
Eigen::Matrix<Derived, size, size> getMatrixAdjugate(const Eigen::Matrix<Derived, size, size>& matrix);
template <typename Derived, int size>
Eigen::Matrix<Derived, size, size> get3x3SymmetricMatrixAdjugate(const Eigen::Matrix<Derived, size, size>& mtx);

template <typename T>
T deg2rad(const T deg);
template <typename T>
T rad2deg(const T rad);
template <typename T>
Eigen::Vector3d latlon2unitVector(const T lat, const T lon);
template <typename T>
Eigen::Vector3d latlon2unitVector(const T crater);
template <typename T>
void operator /(Eigen::Vector3d& vec, T divisor);
template <typename T>
void operator *(Eigen::Vector3d& vec, T scalar);
template <typename T>
Eigen::Vector3d llarToWorld(const T crater, double alt, double rad=1.0);
template <typename R, typename T>
Eigen::Vector3d LLHtoECEF(const R crater, const T alt);
// double vectorNorm(const Eigen::Vector3d& pt);
// void normalizeVector(Eigen::Vector3d& pt);
// template<typename Iter_T>
// long double vectorNorm(Iter_T first, Iter_T last);
// template<typename T>
// T vectorNorm(std::vector<T> vec);
double vdot(const Eigen::Vector3d& pt1, const Eigen::Vector3d& pt2);
double angularPseudodistance(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2);
double angularDistance(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2);
template <typename T>
double latlon_dist(const T crater1, const T crater2);
template <typename T>
double latlon_dist(const T lat1, const T lon1, const T lat2, const T lon2);
template <typename T>
void makeUnique(T& vec);
template <typename T>
std::vector<uint> getRange(std::vector<T> vec);
template <typename T, int size>
bool normalizeDeterminant(Eigen::Matrix<T, size, size>& mtx);
Eigen::Matrix3d crossMatrix(const Eigen::Vector3d&);
void normalizeVector(Eigen::Vector3d&);
void normalizeVector(const Eigen::Vector3d&, Eigen::Vector3d&);
Eigen::Vector3d normalizeVector(const Eigen::Vector3d&);
template <typename T>
T vectorNorm(const Eigen::Vector3d&);

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
Eigen::Vector3d latlon2unitVector(const T lat, const T lon) {
    double lat_rad = deg2rad(lat);
    double lon_rad = deg2rad(lon);

    Eigen::Vector3d point;
    point << cos(lat_rad) * cos(lon_rad),
             cos(lat_rad) * sin(lon_rad),
             sin(lat_rad);

    return point;
}

template <typename T>
Eigen::Vector3d latlon2unitVector(const T crater) {
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

template <typename T>
T vectorNorm(const Eigen::Vector3d& vec) {
  return vec/vec.norm();
}

#endif