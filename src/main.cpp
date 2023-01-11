#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include <map>
#include "crater-id.h"
 
// using namespace std;

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

std::vector<std::tuple<uint, uint, uint>> formTriad(const std::vector<std::tuple<uint, uint>> pairs) {
    std::vector<std::tuple<uint, uint, uint>> combos;
    uint ki, kj; 
    uint idx, iidx, jidx = 0;
    std::vector<std::tuple<uint, uint>>::const_iterator ii, iti, itj = pairs.begin();

    // std::vector<uint>::iterator iti, itj;
    idx = 0;
    for(const auto& [i, j] : pairs) {
        iidx=0;
        do {
            jidx = 0;
            iti = std::find_if(pairs.begin()+idx+iidx, pairs.end(), [&i](const std::tuple<uint,uint>& e) {return std::get<0>(e) == i;});
            do {
                itj = std::find_if(pairs.begin()+idx+iidx+jidx, pairs.end(), [&j](const std::tuple<uint,uint>& e) {return std::get<0>(e) == j;});
                ki = std::get<1>(*iti);
                kj = std::get<1>(*itj);
                std::cout << "--- I:" << i << ", J:" << j << " Ki:" << ki << " Kj:" << kj << std::endl;
                if(    ki == kj 
                    && iti != pairs.end() 
                    && std::get<0>(*iti) == i 
                    && itj != pairs.end() 
                    && std::get<0>(*itj) == j) {
                    // std::tuple<uint, uint, uint> combo = {i, j, ki};
                    combos.push_back({i, j, ki});
                    std::cout << "*** I:" << i << ", J:" << j << " K:" << ki << std::endl;
                }
                jidx++;
            } while(itj != pairs.end() && std::get<0>(*itj) == j);
            iidx++;
        } while(iti != pairs.end() && std::get<0>(*iti) == i);

        idx++;
        ii++;
    }
    return combos;
}
 
int main() {
    std::string fname;
    const std::string degrees = "°";
    // cout<<"Enter crater file name: ";
    // cin>>fname;
    // fname = "/home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv";
    fname = "/home/dqr0509/data/craters/lunar_craters.csv";

    std::vector<lunar_crater> craters = runCraterReader(fname);
    // std::vector<int> choices = {1, 2, 9, 4, 5, 7, 10};
    // std::vector<int> myvector = {10,20,20,20,30,30,20,20,10};

    // printVector(myvector);
    // makeUnique(myvector);
    // printVector(myvector);

    // uint r = 3;
    // std::vector<std::vector<int>> combinations = Combination(choices, r);
    // // printVectorOfVectors(combinations);
    // std::vector<int> rng = getRange(combinations);
    // // printVector(getRange(combinations));

    // // Create a map of strings to integers
    // std::map<int, lunar_crater> mmap;
    // for(size_t idx = 0; idx < craters.size(); idx++) {
    //     mmap[idx] = craters[idx];
    // }
    // std::map<int, lunar_crater>::iterator it = mmap.begin();
    
    // // Iterate through the map and print the elements
    // while (it != mmap.end())
    // {
    //     std::cout << "Key: " << it->first << ", Value: " << it->second << std::endl;
    //     ++it;
    // }
    // std::vector<std::tuple<lunar_crater, lunar_crater>> valids = specialCombination(craters, 10.);
    std::vector<std::tuple<uint, uint>> valids;
    specialCombination(craters, 10., valids);
    
    std::vector<std::tuple<uint, uint, uint>> triads = formTriad(valids);
    uint idx = 0;
    for(const auto& [i, j, k]: triads) {
        std::cout << "IDX: " << idx++ << " | " << i << ", " << j << ", " << k << std::endl;
    }

    // uint count = 0;
    // for(auto& crater : valids) {
    //     std::cout << "-" << count++ << "->"
    //               << " " << std::get<0>(crater)
    //               << "," << std::get<1>(crater)
    //               << std::endl;
    // }
    return 0;
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

std::vector<lunar_crater> runCraterReader(const std::string fname,
                                          const char sep,
                                          const float max_ell,
                                          const float min_arc,
                                          const uint max_n_craters) {
    // martian_crater db_entry;
    lunar_crater   db_entry;
    std::vector<lunar_crater> entries;
    std::string line;

    std::ifstream file (fname, std::ios::in);
    if(file.is_open()) {
        uint count = 0;

        while(getline(file, line)) {
            std::stringstream str(line);
            if(count==0) {
                count++;
                // readLunarCraterEntry(line, sep);
                continue;
            }
            if(readLunarCraterEntry(line, db_entry, sep, max_ell, min_arc)) {
                count++;
                entries.push_back(db_entry);
            }
            if(count > max_n_craters) {
                break;
            }
        }
        file.close();
    }
    else
        std::cout<<"Could not open the file: "<<fname<<std::endl;

    std::cout<<"Database contains "
             <<entries.size()
             <<" entries"
             <<std::endl<<std::endl;
    
    // uint count = 0;
    // uint n_print = 6;
    // for (const auto& crater : entries) {
    //     printLunarCratersInfo(crater);
    //     std::cout<<std::endl;
    //     if(++count >= n_print)
    //         break;
    // }
    return entries;
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

template <typename T, typename N>
std::vector<std::tuple<T, T>> specialCombination(const std::vector<T> choices, 
                                                 const N max_angle_deg) {
    const uint n = choices.size();
    // const float max_angle_deg = 30.;
    const N max_angle_rad = deg2rad(max_angle_deg);
    // trig function are expensive in loops, so use the max dot product
    // angle = acos(dot(a,b)) == cos(angle) = dot(a,b)
    const N min_dot_prod = cos(max_angle_rad);

    Point pt1, pt2;
    std::vector<std::tuple<T, T>> valid_craters;
    std::string lat1, lat2, lon1, lon2;
    T current_choice, next_choice;

    // uint count = 0;
    for(size_t i=0; i<n-1; i++) {
        current_choice = choices[i];
        pt1 = latlon2unitVector(current_choice);
        for(size_t j=i+1; j<n; j++) {
            next_choice = choices[j];
            pt2 = latlon2unitVector(next_choice);
            // if(angularDistance(pt1, pt2) < max_angle_rad) {
            if(angularPseudodistance(pt1, pt2) > min_dot_prod) {
                // count++;
                valid_craters.push_back({current_choice, next_choice});
            }
        }
    }
    // std::cout << "Total number of valid combinations: " << count << std::endl;
    return valid_craters;
}

template <typename T, typename N>
void specialCombination(const std::vector<T> choices, 
                        const N max_angle_deg,
                        std::vector<std::tuple<uint, uint>>& valid_craters) {
    const uint n = choices.size();
    // const float max_angle_deg = 30.;
    const N max_angle_rad = deg2rad(max_angle_deg);
    // trig function are expensive in loops, so use the max dot product
    // angle = acos(dot(a,b)) == cos(angle) = dot(a,b)
    const N min_dot_prod = cos(max_angle_rad);

    Point pt1, pt2;
    // std::vector<std::tuple<T, T>> valid_craters;
    std::string lat1, lat2, lon1, lon2;
    T current_choice, next_choice;

    // uint count = 0;
    for(size_t i=0; i<n-1; i++) {
        current_choice = choices[i];
        pt1 = latlon2unitVector(current_choice);
        for(size_t j=i+1; j<n; j++) {
            next_choice = choices[j];
            pt2 = latlon2unitVector(next_choice);
            // if(angularDistance(pt1, pt2) < max_angle_rad) {
            if(angularPseudodistance(pt1, pt2) > min_dot_prod) {
                // count++;
                valid_craters.push_back({i, j});
            }
        }
    }
    // std::cout << "Total number of valid combinations: " << count << std::endl;
    // return valid_craters;
}

template <typename T>
std::vector<std::vector<T>> Combination(const std::vector<T> choices, const uint r) {
    // https://stackoverflow.com/questions/9430568/generating-combinations-in-c
    const uint n = choices.size();
    std::vector<bool> v(n);
    std::vector<std::vector<T>> output;
    std::vector<T> subvector;
    std::fill(v.begin(), v.begin() + r, true);

    do {
        subvector.clear();
        uint idx = 0;
        for(auto& choice : choices) {
        // for (int i = 0; i < n; ++i) {
            if (v[idx]) {
                subvector.push_back(choice);
            }
            idx++;
        }
        output.push_back(subvector);
    } while (std::prev_permutation(v.begin(), v.end()));
    return output;
}

template <typename T>
void Permutation(std::vector<T> v)
{
    std::sort(v.begin(), v.end());
    do {
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " "));
        std::cout << std::endl;
    } while (std::next_permutation(v.begin(), v.end()));
}

std::ostream& operator<<(std::ostream& os, const lunar_crater& crater) {
    // return os << "Crater: " << crater.crater_id;
    std::string latlon = stringify_latlon(crater);
    return os 
            // << std::fixed << std::setw(10) 
            // << std::setprecision(3) 
            // << std::setfill('_')
            << "ID: " << crater.crater_id
            << "  (" << stringify_latlon(crater) << ")  ";
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

template <typename T>
Point llarToWorld(const T crater, float alt, float rad=1.0) {

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

    // const d = R * c; // in metres
    // std::cout << "A: " << lat1 << "/" << lon1 << "\t";
    // std::cout << "Unit vector: " << point1 << std::endl;
    // std::cout << "B: " << lat2 << "/" << lon2 << "\t";
    // std::cout << "Unit vector: " << point2 << std::endl;
    return angularDistance(point1, point2);
}
