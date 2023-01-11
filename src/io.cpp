#include "io.h"
#include "structs.h"

template <typename T>
void printVector(std::vector<T> vec) {
    // std::cout << "Printing vector: " << std::endl;
    for(auto& idx : vec) {
        std::cout << idx << ", ";
    }
    std::cout << std::endl;
}

template <typename T>
void printVectorOfVectors(std::vector<std::vector<T>> vec) {
    std::cout << "Printing vector of vectors: " << std::endl;
    for(auto& combo : vec) {
        printVector(combo);
    }
    std::cout << std::endl;
}

bool readLunarCraterEntry(const std::string entry, 
                          lunar_crater& crater,
                          const char sep, 
                          const float max_ell,
                          const float min_arc,
                          const float min_diam,
                          const float max_diam) {
    // lunar_crater crater;
    std::string lat, lon, diam, ecc, n_pts;
    std::stringstream str(entry);
    std::string token;
    uint count = 0;
    while(getline(str, token, sep)) {
        if(token.empty())
            return false;
        switch(count) {
            case CRATER_ID:
                crater.crater_id = token;
                break;
            case LAT_CIRC_IMG:
                crater.lat = stof(token);
                break;
            case LON_CIRC_IMG:
                crater.lon = stof(token);
                break;
            case DIAM_CIRC_IMG:
                crater.diam = stof(token);
                break;
            case DIAM_ELLI_ECCEN_IMG:
                crater.ecc = stof(token);
                break;
            case DIAM_ELLI_ELLIP_IMG:
                crater.ell = stof(token);
                break;
            case ARC_IMG:
                crater.arc = stof(token);
                break;
            default:
                break;
        }
        count++;
    }
    return crater.ell < max_ell 
        && crater.arc > min_arc
        && crater.diam < max_diam
        && crater.diam > min_diam;
}

bool readMartianCraterEntry(const std::string entry) {
    return false;
}

template <typename T>
std::tuple<float, char> stringify_lat(const T crater) {
    return stringify_lat(crater.lat);
}

template <typename T>
std::tuple<float, char> stringify_lon(const T crater) {
    return stringify_lon(crater.lon);
}

std::tuple<float, char> stringify_lat(const float lat) {
    char n_or_s;
    float re_lat;
    // tuple<float, std::string> re_lat;
    n_or_s = (lat > 0)?'N':'S';
    re_lat = abs(lat);
    return {re_lat, n_or_s};
}

std::tuple<float, char> stringify_lon(const float lon) {
    char e_or_w;
    float re_lon = lon;
    if(re_lon > 180)
        re_lon -= 360;
    e_or_w = (re_lon < 0)?'W':'E';
    // e_or_w = (lon < 0 || lon > 180)?'W':'E';        
    return {abs(re_lon), e_or_w};
}

std::string stringify_latlon(const float lat, const float lon) {
    float re_lat, re_lon;
    char n_or_s, e_or_w;
    std::string str_lat, str_lon;
    std::tie(re_lat, n_or_s) = stringify_lat(lat);
    std::tie(re_lon, e_or_w) = stringify_lon(lon);
    str_lat = std::to_string(re_lat) + "°" + n_or_s;
    str_lon = std::to_string(re_lon) + "°" + e_or_w;    
    return (str_lat + " " + str_lon);
}

template <typename T>
std::string stringify_latlon(const T crater) {
    return stringify_latlon(crater.lat, crater.lon);
}

void printLunarCratersInfo(const lunar_crater crater) {
    // const char cout_sep = '\n';
    const std::string cout_sep = "\n\t";

    std::string latlon = stringify_latlon(crater);
    std::cout 
         << "CRATER_ID: " << crater.crater_id << cout_sep
         << "LAT_CIRC_IMG/LON_CIRC_IMG: " << latlon << cout_sep 
         << "DIAM_CIRC_IMG:       " << crater.diam << "km" << cout_sep 
         << "DIAM_ELLI_ECCEN_IMG: " << crater.ecc << cout_sep 
         << "DIAM_ELLI_ELLIP_IMG: " << crater.ell << cout_sep 
         << "ARC_IMG:             " << crater.arc << cout_sep
         << std::endl;
}

bool readLunarCraterEntry(const std::string entry, const char sep) {
    // lunar_crater crater;
    std::string lat, lon, diam, ecc, n_pts;
    std::stringstream str(entry);
    std::string token;
    uint count = 0;
    while(getline(str, token, sep)) {
        if(token.empty())
            return false;
        switch(count++) {
            case CRATER_ID:
            case LAT_CIRC_IMG:
            case LON_CIRC_IMG:
            case DIAM_CIRC_IMG:
            case DIAM_ELLI_ECCEN_IMG:
            case DIAM_ELLI_ELLIP_IMG:
            case ARC_IMG:
                std::cout<<token<<sep;
            default:
                break;
        }
    }
    std::cout<<std::endl;
    return true;
}