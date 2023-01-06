#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
// #include "include/crater-id.h"
 
using namespace std;

struct martian_crater {
    string CRATER_ID;
    float  LATITUDE_CIRCLE_IMAGE;
    float  LONGITUDE_CIRCLE_IMAGE;
    float  DIAM_CIRCLE_IMAGE;
    float  DIAM_ELLIPSE_ECCEN_IMAGE;
    uint    PTS_USED_RIM_IMAGE;
    string CRATER_NAME;
};

struct lunar_crater {
    string CRATER_ID;
	float  LAT_CIRC_IMG;
	float  LON_CIRC_IMG;
	float  DIAM_CIRC_IMG;
	float  DIAM_CIRC_SD_IMG;
	float  DIAM_ELLI_ECCEN_IMG;
	float  DIAM_ELLI_ELLIP_IMG;
	float  ARC_IMG;
	uint   PTS_RIM_IMG;
};

bool readMartianCraterEntry(const string entry, const char sep) {
    martian_crater crater;
    string lat, lon, diam, ecc, n_pts;
    stringstream str(entry);
    getline(str, crater.CRATER_ID, sep);
    getline(str, lat, sep);
    getline(str, lon, sep);
    getline(str, diam, sep);
    getline(str, ecc,  sep);
    getline(str, n_pts,sep);
    getline(str, crater.CRATER_NAME, sep);
    try {
        crater.LATITUDE_CIRCLE_IMAGE = stof(lat);
        crater.LONGITUDE_CIRCLE_IMAGE = stof(lon);
        crater.DIAM_CIRCLE_IMAGE = stof(diam);
        crater.DIAM_ELLIPSE_ECCEN_IMAGE = stof(ecc);
        crater.PTS_USED_RIM_IMAGE = stof(n_pts);
    } catch (const exception& e) {return false;}
    return crater.DIAM_ELLIPSE_ECCEN_IMAGE < 0.3;
}

bool readLunarCraterEntry(const string entry) {
    return false;
}
 
int main() {
    string fname;
    // cout<<"Enter crater file name: ";
    // cin>>fname;
    fname = "/home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv";
    
    martian_crater db_entry;
    vector<vector<string>> content;
    vector<string> row;
    vector<martian_crater> entries;
    string line, word;
    // string crater_id, crater_name;
    string lat, lon, diam, ecc, n_pts;

    // /home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv
    char sep = ',';
    // char cout_sep = '\t';
    
    ifstream file (fname, ios::in);
    if(file.is_open()) {
        uint count = 0;
        while(getline(file, line)) {
            row.clear();
            count++;
            stringstream str(line);
            if(count<2) {
                cout<<line<<endl;
                while(getline(str, word, sep))
                    row.push_back(word);
                continue;
            }
            getline(str, db_entry.CRATER_ID, sep);
            getline(str, lat, sep);
            getline(str, lon, sep);
            getline(str, diam, sep);
            getline(str, ecc,  sep);
            getline(str, n_pts,sep);
            getline(str, db_entry.CRATER_NAME, sep);
            try {
            db_entry.LATITUDE_CIRCLE_IMAGE = stof(lat);
            db_entry.LONGITUDE_CIRCLE_IMAGE = stof(lon);
            db_entry.DIAM_CIRCLE_IMAGE = stof(diam);
            db_entry.DIAM_ELLIPSE_ECCEN_IMAGE = stof(ecc);
            db_entry.PTS_USED_RIM_IMAGE = stof(n_pts);
            } catch (const exception& e) {continue;}
            if(db_entry.DIAM_ELLIPSE_ECCEN_IMAGE > 0.9)
                entries.push_back(db_entry);
            // cout << db_entry.CRATER_ID << cout_sep 
            //      << db_entry.LATITUDE_CIRCLE_IMAGE << cout_sep 
            //      << db_entry.LONGITUDE_CIRCLE_IMAGE << cout_sep 
            //      << db_entry.DIAM_CIRCLE_IMAGE << cout_sep 
            //      << db_entry.DIAM_ELLIPSE_ECCEN_IMAGE << cout_sep 
            //      << db_entry.PTS_USED_RIM_IMAGE << cout_sep 
            //      << db_entry.CRATER_NAME << endl;
        }
    }
    else
        cout<<"Could not open the file: "<<fname<<endl;
    uint count = 0;
    for (const auto& cline : content) {
        for (const auto& data : cline) {
            cout<<data<<"  \t";
        }
        cout<<endl;
        if(count++> 10)
            break;
    }
    
    return 0;
}