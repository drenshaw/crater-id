#include "vector_math.h"

template <typename T>
Eigen::Vector3f llarToWorld(const T crater, float alt, float rad) {

    const T lat = crater.lat;
    const T lon = crater.lon;
    // see: http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html
    T f  = 0;                              // flattening
    T ls = atan2(pow(1 - f, 2) * tan(lat));    // lambda

    Eigen::Vector3f point;
    point << rad * cos(ls) * cos(lon) + alt * cos(lat) * cos(lon),
             rad * cos(ls) * sin(lon) + alt * cos(lat) * sin(lon),
             rad * sin(ls) + alt * sin(lat);

    return point;
}

template <typename R, typename T>
Eigen::Vector3f LLHtoECEF(const R crater, const T alt) {
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

    Eigen::Vector3f point;
    point << (rad * C + alt)*cosLat * cos(lon),
             (rad * C + alt)*cosLat * sin(lon),
             (rad * S + alt)*sinLat;

    return point;
}

float vdot(const Eigen::Vector3f& point1, const Eigen::Vector3f& point2) {
    return point1.dot(point2);
}

float angularPseudodistance(const Eigen::Vector3f& point1, const Eigen::Vector3f& point2) {
    return vdot(point1, point2);
}

float angularDistance(const Eigen::Vector3f& point1, const Eigen::Vector3f& point2) {
    return acos(vdot(point1, point2));
}

template <typename T>
float latlon_dist(const T crater1, const T crater2) {
    return latlon_dist(crater1.lat, crater1.lon, crater2.lat, crater2.lon);
}

template <typename T>
float latlon_dist(const T lat1, const T lon1, const T lat2, const T lon2) {
    Eigen::Vector3f point1 = latlon2unitVector(lat1, lon1);
    // normalizeVector(point1);
    Eigen::Vector3f point2 = latlon2unitVector(lat2, lon2);
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

Eigen::Matrix3f normalizeDeterminant(const Eigen::Matrix3f& mtx) {
    // using approach from Matrix Cookbook
    // https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf
    // get the determinant
    uint ncol = mtx.cols();
    float det_mtx = mtx.determinant();
    if(abs(det_mtx)<1e-20) {
        std::cerr << "Matrix is singular/nearly singular." << std::endl;
        return mtx;
    }
    // we want the determinant of A to be 1
    // 1 = det(c*A) = d^n*det(A), where n=3 for 3x3 matrix
    // so solve the equation d^3 - 1/det(A) = 0
    // d = roots([1 0 0 -1/detC]);
    // strip out imaginary roots
    // d = d(imag(d)==0);
    // assert that we only have one real root
    // assert(length(d)==1)
    float d = pow(abs(1./det_mtx), 1./ncol);
    // d*mtx should now have a determinant of 1
    Eigen::Matrix3f mtx_normalized = d*mtx;
    int sign = signbit(mtx_normalized.determinant());
    // assert(abs(det(A_norm)-1)<sqrt(eps))
    return sign? -mtx_normalized : mtx_normalized;
}

Eigen::Matrix3f crossMatrix(const Eigen::Vector3f& v) {
    Eigen::Matrix3f v_cross(3, 3);
    v_cross << 0, -v(2), v(1),
              v(2), 0, -v(0),
              -v(1), v(0), 0;
    return v_cross;
}