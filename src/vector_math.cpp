#include "vector_math.h"

template <typename T>
Point llarToWorld(const T crater, float alt, float rad) {

    const T lat = crater.lat;
    const T lon = crater.lon;
    // see: http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html
    T f  = 0;                              // flattening
    T ls = atan2(pow(1 - f, 2) * tan(lat));    // lambda

    Point point;
    point.x = rad * cos(ls) * cos(lon) + alt * cos(lat) * cos(lon);
    point.y = rad * cos(ls) * sin(lon) + alt * cos(lat) * sin(lon);
    point.z = rad * sin(ls) + alt * sin(lat);

    return point;
}

template <typename R, typename T>
Point LLHtoECEF(const R crater, const T alt) {
    // see http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html

    const T lat = crater.lat;
    const T lon = crater.lon;
    const T rad = 6378137.0;        // Radius of the Earth (in meters)
    const T f = 1.0/298.257223563;  // Flattening factor WGS84 Model
    T cosLat = cos(lat);
    T sinLat = sin(lat);
    T FF     = pow(1.0-f, 2);
    T C      = 1/sqrt(pow(cosLat,2) + FF * pow(sinLat,2));
    T S      = C * FF;

    Point point;
    point.x = (rad * C + alt)*cosLat * cos(lon);
    point.y = (rad * C + alt)*cosLat * sin(lon);
    point.z = (rad * S + alt)*sinLat;

    return point;
}

float vectorNorm(const Point pt) {
    return sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z);
}

void normalizeVector(Point& pt) {
    float norm = vectorNorm(pt);
    pt.x /= norm;
    pt.y /= norm;
    pt.z /= norm;
}

// template<typename Iter_T>
// long double vectorNorm(Iter_T first, Iter_T last) {
//   return sqrt(inner_product(first, last, first, 0.0L));
// }

// template<typename T>
// void vectorNorm(std::vector<T> vec) {
//     // std::cout << "Norm: " << vectorNorm(vec.begin(), vec.end());
// //   return sqrt(inner_product(vec.begin(), vec.end(), vec.begin(), 0.0L));
//     // return 0.;
// }

float dot(const Point pt1, const Point pt2) {
    return pt1.x*pt2.x + pt1.y*pt2.y + pt1.z*pt2.z;
}

float angularPseudodistance(const Point point1, const Point point2) {
    return dot(point1, point2);
}

float angularDistance(const Point point1, const Point point2) {
    return acos(dot(point1, point2));
}

template <typename T>
float latlon_dist(const T crater1, const T crater2) {
    return latlon_dist(crater1.lat, crater1.lon, crater2.lat, crater2.lon);
}

template <typename T>
float latlon_dist(const T lat1, const T lon1, const T lat2, const T lon2) {
    Point point1 = latlon2unitVector(lat1, lon1);
    // normalizeVector(point1);
    Point point2 = latlon2unitVector(lat2, lon2);
    // normalizeVector(point2);
    return angularDistance(point1, point2);
}

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