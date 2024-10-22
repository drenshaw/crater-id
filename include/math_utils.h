#pragma once
// TODO: rename this to a better name like `utils` or something

#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <eigen3/Eigen/Dense>
// #include <boost/log/core.hpp>
// #include <boost/log/expressions.hpp>
// #include <boost/log/trivial.hpp>
// #include "crater-id.h"

#define R_MOON 1737.4
const double EPS = (10 * std::numeric_limits<double>::epsilon());

double getCofactor(const Eigen::MatrixXd& matrix, int p, int q);
Eigen::MatrixXd cofactor(const Eigen::MatrixXd& matrix);
Eigen::MatrixXd adjugate(const Eigen::MatrixXd&);
Eigen::Matrix3d symmetricAdjugate(const Eigen::Matrix3d&);

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
// double vectorNorm(const Eigen::Vector3d& pt);
// void normalizeVector(Eigen::Vector3d& pt);
// template<typename Iter_T>
// long double vectorNorm(Iter_T first, Iter_T last);
// template<typename T>
// T vectorNorm(std::vector<T> vec);
double vdot(const Eigen::Vector3d& pt1, const Eigen::Vector3d& pt2);
double getPseudoAngleBetweenVectors(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2);
double getAngleBetweenVectors(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2);
Eigen::Vector3d getAxisNormalToVectors(const Eigen::Vector3d& vec1, const Eigen::Vector3d& vec2);

Eigen::Matrix3d crossMatrix(const Eigen::Vector3d&);
template <typename T>
T vectorNorm(const Eigen::Vector3d&);
Eigen::Matrix3d getENUFrame(const Eigen::Vector3d&);
Eigen::Matrix3d getENUFrame(const double, const double);
Eigen::Vector3d getNorthPoleUnitVector();
void GetNorthPoleUnitVector(Eigen::Vector3d&);
Eigen::Quaterniond eulerToQuaternion(const double roll, const double pitch, const double yaw);
void eulerToDCM(const double roll,
                const double pitch,
                const double yaw,
                Eigen::Matrix3d& dcm);    

bool vectorContainsNaN(const Eigen::Vector3d& eV);
void convertEigenVectorToVector(const Eigen::Vector3d& eig, std::vector<double>& vec);


/**** Template definitions ****/

template <typename Derived>
bool normalizeDeterminant(Eigen::MatrixBase<Derived>& mtx) {
  // using approach from Matrix Cookbook
  // https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf
  // get the determinant
  if(mtx.isApprox(Eigen::MatrixBase<Derived>::Zero())) {
    std::cerr << __func__ << ": Matrix is zero.\n";
    throw std::runtime_error("Matrix is zero");
  }
  uint ncol = mtx.cols();
  assert(mtx.cols() == mtx.rows());
  // uint nrow = mtx.rows();
  // T mtx_normalized(nrow, ncol);
  //get location of maximum
  // Eigen::Index maxRow, maxCol;
  double maxVal = 1/mtx.maxCoeff();
  if(std::isnan(maxVal)) {
    // BOOST_LOG_TRIVIAL(warning) << "Matrix is singular.";
    std::cerr << __func__ << "Max coefficient is 0.\n";
    return false;
  }
  mtx *= maxVal;
  double det_mtx = mtx.determinant();
  if(det_mtx == 0) {
    // BOOST_LOG_TRIVIAL(warning) << "Matrix is singular.";
    std::cerr << "Matrix is singular/nearly singular." << std::endl;
    mtx *= maxVal;
    return false;
  }
  // we want the determinant of A to be 1
  // 1 = det(c*A) = d^n*det(A), where n=3 for 3x3 matrix
  // so solve the equation d^3 - 1/det(A) = 0
  double d = pow(abs(1./det_mtx), 1./ncol);
  // d*mtx should now have a determinant of +1
  mtx *= d;
  bool sign = std::signbit(mtx.determinant()); // true if negative
  // assert(abs(det(A_norm)-1)<sqrt(eps))
  mtx *= sign ? -1.0 : 1.0;
  return true;
}

template <typename T, size_t SIZE>
void convertEigenVectorToArray(const Eigen::Vector3d& eig, std::array<T, SIZE>& arr){
  Eigen::Vector3d::Map(&arr[0], eig.size()) = eig;
}

template <typename T, size_t SIZE>
bool vectorContainsNaN(const std::array<T, SIZE>& vec) {
  return std::any_of(vec.begin(), vec.end(), [](T i){return std::isnan(i);});
}
template <typename T, size_t SIZE>
bool vectorContainsNaN(const std::vector<T>& vec) {
  return std::any_of(vec.begin(), vec.end(), [](T i){return std::isnan(i);});
}
template <typename T, size_t SIZE>
bool vectorContainsNaN(const Eigen::Vector3d& eV) {
  std::array<T, SIZE> vec;
  convertEigenVectorToArray(eV, vec);
  return vectorContainsNaN(vec);
}

template <typename T, typename A>
int arg_max(std::vector<T, A> const& vec) {
  return static_cast<int>(
    std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
}

template <typename T, typename A>
int arg_min(std::vector<T, A> const& vec) {
  return static_cast<int>(
    std::distance(vec.begin(), min_element(vec.begin(), vec.end())));
}

template <typename T, size_t SIZE>
int arg_max(std::array<T, SIZE> const& arr) {
  return static_cast<int>(
    std::distance(arr.begin(), std::max_element(arr.begin(), arr.end())));
}

template <typename T, size_t SIZE>
int arg_min(std::array<T, SIZE> const& arr) {
  return static_cast<int>(
    std::distance(arr.begin(), std::min_element(arr.begin(), arr.end())));
}

template <typename T>
double deg2rad(const T deg) {
  return static_cast<double>(deg) * M_PI/180.;
}
template <typename T>
double rad2deg(const T rad) {
  return static_cast<double>(rad) * 180. / M_PI;
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


template <typename T>
double getAngleBetweenLatLon(const T crater1, const T crater2) {
  return getAngleBetweenLatLon(crater1.lat, crater1.lon, crater2.lat, crater2.lon);
}

template <typename T>
double getAngleBetweenLatLon(const T lat1, const T lon1, const T lat2, const T lon2) {
  Eigen::Vector3d point1 = latlon2bearing(lat1, lon1);
  // normalizeVector(point1);
  Eigen::Vector3d point2 = latlon2bearing(lat2, lon2);
  // normalizeVector(point2);
  return getAngleBetweenVectors(point1, point2);
}

double calculateCraterRimFromRadius(const double radius);
Eigen::Vector3d latlonrad2CraterRim(const double lat, const double lon, const double radius);
Eigen::Vector3d latlonalt(const double lat, const double lon, const double altitude);

template <typename T>
void makeUnique(T& vec) {
  std::sort(vec.begin(), vec.end());
  std::vector<int>::iterator it;
  it = std::unique(vec.begin(), vec.end());
  vec.resize(std::distance(vec.begin(), it));
}

template <typename T>
std::vector<uint> getRange(std::vector<T> vec) {
  // std::vector<uint> vec_range(vec.size());
  std::vector<uint> vec_range;
  vec_range.reserve(vec.size());
  for(size_t idx = 0; idx < vec.size(); idx++) {
    // vec_range[idx] = idx;
    vec_range.push_back(idx);
  }
  return vec_range;
}

// TODO: Use LLHtoECEF instead for the time being for consistency
// TODO: Sometimes, this and LLHtoECEF get negatives of each other
template <typename T>
Eigen::Vector3d llarToECEF(const T lat, const T lon, const double alt, const double radius=1., const double flattening=0.) {
  // see: http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html
  // TODO: ensure that flattening value is valid
  assert(flattening > 0.0 && flattening < 1.0);
  T ls = atan(pow(1 - flattening, 2.0) * tan(lat));  // lambda

  Eigen::Vector3d point;
  point <<  radius * cos(ls) * cos(lon) + alt * cos(lat) * cos(lon),
            radius * cos(ls) * sin(lon) + alt * cos(lat) * sin(lon),
            radius * sin(ls) + alt * sin(lat);

  return point;
}

// TODO: Use LLHtoECEF instead for the time being for consistency
template <typename T>
Eigen::Vector3d llarToECEF(const T crater, const double alt, const double radius=1., const double flattening=0.) {
  const T lat = crater.lat;
  const T lon = crater.lon;
  // see: http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html
  return llarToECEF(lat, lon, alt, radius, flattening);
}

template <typename T>
Eigen::Vector3d LLHtoECEF(const T lat, const T lon, const T alt) {
    // see http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html
  const T radius = 6378137.0;     // Radius of the Earth (in meters)
  const T f = 1.0/298.257223563;  // Flattening factor WGS84 Model
  T cosLat = cos(lat);
  T sinLat = sin(lat);
  T FF     = pow(1.0-f, 2);
  T C      = 1/sqrt(pow(cosLat,2.0) + FF * pow(sinLat,2));
  T S      = C * FF;

  Eigen::Vector3d point;
  point <<  (radius * C + alt)*cosLat * cos(lon),
            (radius * C + alt)*cosLat * sin(lon),
            (radius * S + alt)*sinLat;

  return point;
}

template <typename R, typename T>
Eigen::Vector3d LLHtoECEF(const R crater, const T alt) {
    // see http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html
  const T lat = crater.lat;
  const T lon = crater.lon;
  return LLHtoECEF(lat, lon, alt);
}
