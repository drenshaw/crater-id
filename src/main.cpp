#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include "crater-id.h"
 
// using namespace std;
 
int main() {
    std::string fname;
    // cout<<"Enter crater file name: ";
    // cin>>fname;
    // fname = "/home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv";
    fname = "/home/dqr0509/data/craters/lunar_craters.csv";

    std::vector<lunar_crater> craters = runCraterReader(fname);
    std::vector<std::vector<lunar_crater>> combinations;
    std::vector<int> choices = {1, 2, 9, 4, 5, 7, 10};
    // std::vector<char> choices = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'};
    // Permutation(choices);


    int r = 3;
    combinations = Combination(craters, r);
    return 0;
}

bool readLunarCraterEntry(  const std::string entry, 
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
    
    uint count = 0;
    uint n_print = 6;
    for (const auto& crater : entries) {
        printLunarCratersInfo(crater);
        std::cout<<std::endl;
        if(++count >= n_print)
            break;
    }
    return entries;
}

bool readMartianCraterEntry(const std::string entry) {
    return false;
}

std::tuple<float, char> reinterpret_lat(const float lat) {
    char n_or_s;
    float re_lat;
    // tuple<float, std::string> re_lat;
    n_or_s = (lat > 0)?'N':'S';
    re_lat = abs(lat);
    return {re_lat, n_or_s};
}

std::tuple<float, char> reinterpret_lon(const float lon) {
    char e_or_w;
    float re_lon;
    e_or_w = (lon < 0 || lon > 180)?'W':'E';
    re_lon = abs(lon);
    return {re_lon, e_or_w};
}

std::tuple<std::string, std::string> reinterpret_latlon(const float lat, 
                                                        const float lon) {
    float re_lat, re_lon;
    char n_or_s, e_or_w;
    std::string str_lat, str_lon;
    std::tie(re_lat, n_or_s) = reinterpret_lat(lat);
    std::tie(re_lon, e_or_w) = reinterpret_lon(lon);
    str_lat = std::to_string(re_lat) + "°" + n_or_s;
    str_lon = std::to_string(re_lon) + "°" + e_or_w;    
    return {str_lat, str_lon};
}

void printLunarCratersInfo(const lunar_crater crater) {
    // const char cout_sep = '\n';
    const std::string cout_sep = "\n\t";
    std::string lat, lon;
    tie(lat, lon) = reinterpret_latlon(crater.lat, crater.lon);
    std::cout << "CRATER_ID: " << crater.crater_id << cout_sep
         << "LAT_CIRC_IMG:        " << lat << cout_sep 
         << "LON_CIRC_IMG:        " << lon << cout_sep 
         << "DIAM_CIRC_IMG:       " << crater.diam << "km" << cout_sep 
         << "DIAM_ELLI_ECCEN_IMG: " << crater.ecc << cout_sep 
         << "DIAM_ELLI_ELLIP_IMG: " << crater.ell << cout_sep 
         << "ARC_IMG:             " << crater.arc << std::endl;
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

template <typename T>
std::vector<std::vector<T>> Combination(const std::vector<T> choices, const int r) {
    const int n = choices.size();
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
    // for(auto& vec : output) {
    //     for(auto& element : vec) {
    //         std::cout << element << " | ";
    //     }
    //     std::cout << std::endl;
    // }
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
    return os << "Crater: " << crater.crater_id;
    // return os << std::fixed << std::setw(10) 
    //           << std::setprecision(0) 
    //           << std::setfill(' ')
    //           << "Crater: " << crater.crater_id
    //           << "\tdiam: " << crater.diam << "km\t";
}
