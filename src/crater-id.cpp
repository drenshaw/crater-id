#include "crater-id.h"

std::ostream& operator<<(std::ostream& os, const lunar_crater& crater) {
    // return os << "Crater: " << crater.crater_id;
    std::string latlon = stringify_latlon(crater);
    return os 
            // << std::fixed << std::setw(10) 
            // << std::setprecision(3) 
            // << std::setfill('_')
            << "ID: " << crater.crater_id
            << "  (" << latlon << ")";
            // << "\tdiam: " << crater.diam << "km ";
}

std::ostream& operator<<(std::ostream& os, const Point& point) {
    return os 
            // << std::fixed << std::setw(10) 
            // << std::setprecision(4) 
            // << std::setfill('_')
            << " ( " << point.x
            << ", " << point.y
            << ", " << point.z << " ) ";
}
 
std::ostream& operator<<(std::ostream& os, const Rect& rectangle) {
    return os 
        // << std::fixed << std::setw(10) 
        // << std::setprecision(3) 
        // << std::setfill('_')
        << "->Rect: [" 
        << rectangle.min[0] << "," 
        << rectangle.min[1] << "] ["
        << rectangle.max[0] << "," 
        << rectangle.max[1] << "] ";
}
 
std::ostream& operator<<(std::ostream& os, const box& bbox) {
    return os 
        // << std::fixed << std::setw(10) 
        // << std::setprecision(3) 
        // << std::setfill('_')
        << "->Box: [" 
        << bg::get<bg::min_corner, 0>(bbox) << "," 
        << bg::get<bg::min_corner, 1>(bbox) << "] ["
        << bg::get<bg::max_corner, 0>(bbox) << "," 
        << bg::get<bg::max_corner, 1>(bbox) << "] ";
}

double calculateCraterRimFromRadius(const double radius) {
    return sqrt(pow(R_MOON, 2) - pow(radius, 2));
}