#pragma once

// #include <iostream>
// #include <math.h>
// #include <boost/geometry.hpp>
// #include <boost/geometry/index/rtree.hpp>
// #include "structs.h"
// #include "io.h"

// namespace bg = boost::geometry;
// namespace bgi = boost::geometry::index;

// std::ostream& operator<<(std::ostream& os, const Point& point);

// typedef bg::model::point<float, NDIM, bg::cs::cartesian> point;
// typedef bg::model::box<point> box;
// typedef std::pair<box, size_t> value;

// struct Rect
// {
//   Rect()  {}

//   Rect(int a_minX, int a_minY, int a_maxX, int a_maxY)
//   {
//     min[0] = a_minX;
//     min[1] = a_minY;

//     max[0] = a_maxX;
//     max[1] = a_maxY;
//   }
//   Rect(int X, int Y)
//   {
//     min[0] = X;
//     min[1] = Y;

//     max[0] = X;
//     max[1] = Y;
//   }

//   int min[2];
//   int max[2];

//   box calculate_bounding_box() {
//     point point1(min[0], min[1]);
//     point point2(max[0], min[1]);
//     box bbox(point1, point2);
//     return bbox;
//   }

// };
