#include "crater-id.h"

std::ostream& operator<<(std::ostream& os, const lunar_crater& crater) {
    // return os << "Crater: " << crater.crater_id;
    std::string latlon = stringify_latlon(crater);
    return os 
            // << std::fixed << std::setw(10) 
            // << std::setprecision(3) 
            // << std::setfill('_')
            << "ID: " << crater.crater_id
            << "  (" << stringify_latlon(crater) << ")";
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