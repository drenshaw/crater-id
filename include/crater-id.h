#ifndef CRATER_ID_H
#define CRATER_ID_H

#include <iostream>
#include <math.h>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "structs.h"
#include "io.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


#define NDIM 2
#define R_MOON 1737.4

std::ostream& operator<<(std::ostream& os, const lunar_crater& crater);
std::ostream& operator<<(std::ostream& os, const Point& point);

double calculateCraterRimFromRadius(const double radius);

typedef bg::model::point<float, NDIM, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<box, size_t> value;

struct Rect
{
  Rect()  {}

  Rect(int a_minX, int a_minY, int a_maxX, int a_maxY)
  {
    min[0] = a_minX;
    min[1] = a_minY;

    max[0] = a_maxX;
    max[1] = a_maxY;
  }
  Rect(int X, int Y)
  {
    min[0] = X;
    min[1] = Y;

    max[0] = X;
    max[1] = Y;
  }

  int min[2];
  int max[2];

  box calculate_bounding_box() {
    point point1(min[0], min[1]);
    point point2(max[0], min[1]);
    box bbox(point1, point2);
    return bbox;
  }

};

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
// Point latlon2bearing(const T lat, const T lon);
// template <typename T>
// Point latlon2bearing(const T crater);
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

#endif