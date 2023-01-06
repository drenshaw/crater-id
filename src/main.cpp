#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include <tuple>
#include "crater-id.h"
 
using namespace std;
 
int main() {
    string fname;
    // cout<<"Enter crater file name: ";
    // cin>>fname;
    // fname = "/home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv";
    fname = "/home/dqr0509/data/craters/lunar_craters.csv";

    // // /home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv

    vector<lunar_crater> craters = runCraterReader(fname);
    
    return 0;
}

bool readLunarCraterEntry(  const string entry, 
                            lunar_crater& crater,
                            const char sep, 
                            const float max_ell,
                            const float min_arc,
                            const float min_diam,
                            const float max_diam) {
    // lunar_crater crater;
    string lat, lon, diam, ecc, n_pts;
    stringstream str(entry);
    string token;
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

vector<lunar_crater> runCraterReader(const string fname,
                                     const char sep,
                                     const float max_ell,
                                     const float min_arc,
                                     const uint max_n_craters) {
    // martian_crater db_entry;
    lunar_crater   db_entry;
    vector<lunar_crater> entries;
    string line;

    ifstream file (fname, ios::in);
    if(file.is_open()) {
        uint count = 0;

        while(getline(file, line)) {
            stringstream str(line);
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
        cout<<"Could not open the file: "<<fname<<endl;

    cout<<"Database contains "<<entries.size()<<" entries"<<endl<<endl;
    
    uint count = 0;
    uint n_print = 6;
    for (const auto& crater : entries) {
        printLunarCratersInfo(crater);
        cout<<endl;
        if(++count >= n_print)
            break;
    }
    return entries;
}

bool readMartianCraterEntry(const string entry) {
    return false;
}

tuple<float, char> reinterpret_lat(const float lat) {
    char n_or_s;
    float re_lat;
    // tuple<float, string> re_lat;
    n_or_s = (lat > 0)?'N':'S';
    re_lat = abs(lat);
    return {re_lat, n_or_s};
}

tuple<float, char> reinterpret_lon(const float lon) {
    char e_or_w;
    float re_lon;
    e_or_w = (lon < 0 || lon > 180)?'W':'E';
    re_lon = abs(lon);
    return {re_lon, e_or_w};
}

tuple<string, string> reinterpret_latlon(const float lat, const float lon) {
    float re_lat, re_lon;
    char n_or_s, e_or_w;
    string str_lat, str_lon;
    tie(re_lat, n_or_s) = reinterpret_lat(lat);
    tie(re_lon, e_or_w) = reinterpret_lon(lon);
    str_lat = to_string(re_lat) + "°" + n_or_s;
    str_lon = to_string(re_lon) + "°" + e_or_w;    
    return {str_lat, str_lon};
}

void printLunarCratersInfo(const lunar_crater crater) {
    // const char cout_sep = '\n';
    const string cout_sep = "\n\t";
    string lat, lon;
    tie(lat, lon) = reinterpret_latlon(crater.lat, crater.lon);
    cout << "CRATER_ID: " << crater.crater_id << cout_sep
         << "LAT_CIRC_IMG:        " << lat << cout_sep 
         << "LON_CIRC_IMG:        " << lon << cout_sep 
         << "DIAM_CIRC_IMG:       " << crater.diam << "km" << cout_sep 
         << "DIAM_ELLI_ECCEN_IMG: " << crater.ecc << cout_sep 
         << "DIAM_ELLI_ELLIP_IMG: " << crater.ell << cout_sep 
         << "ARC_IMG:             " << crater.arc << endl;
}

bool readLunarCraterEntry(const string entry, const char sep) {
    // lunar_crater crater;
    string lat, lon, diam, ecc, n_pts;
    stringstream str(entry);
    string token;
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
                cout<<token<<sep;
            default:
                break;
        }
    }
    cout<<endl;
    return true;
}
