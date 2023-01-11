#pragma once

#include <math.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include "crater-id.h"

template <typename T>
T deg2rad(const T deg);
template <typename T>
T rad2deg(const T rad);
template <typename T>
Point latlon2unitVector(const T lat, const T lon);
template <typename T>
Point latlon2unitVector(const T crater);
template <typename T>
Point llarToWorld(const T crater, float alt, float rad=1.0);
template <typename R, typename T>
Point LLHtoECEF(const R crater, const T alt);
float vectorNorm(const Point pt);
void normalizeVector(Point& pt);
template<typename Iter_T>
long double vectorNorm(Iter_T first, Iter_T last);
template<typename T>
T vectorNorm(std::vector<T> vec);
float dot(const Point pt1, const Point pt2);
float angularPseudodistance(const Point point1, const Point point2);
float angularDistance(const Point point1, const Point point2);
template <typename T>
float latlon_dist(const T crater1, const T crater2);
template <typename T>
float latlon_dist(const T lat1, const T lon1, const T lat2, const T lon2);
template <typename T>
void makeUnique(T& vec);
template <typename T>
std::vector<uint> getRange(std::vector<T> vec);

/**** Template definitions ****/

template <typename T>
T deg2rad(const T deg) {
    return deg * M_PI/180.;
}
template <typename T>
T rad2deg(const T rad) {
    return rad * 180. / M_PI;
}

template <typename T>
Point latlon2unitVector(const T lat, const T lon) {
    float lat_rad = deg2rad(lat);
    float lon_rad = deg2rad(lon);

    Point point;
    point.x = cos(lat_rad) * cos(lon_rad);
    point.y = cos(lat_rad) * sin(lon_rad);
    point.z = sin(lat_rad);

    return point;
}

template <typename T>
Point latlon2unitVector(const T crater) {
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