#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <tuple>
#include "structs.h"

template <typename T>
void printVector(std::vector<T>);
template <typename T>
void printVectorOfVectors(std::vector<std::vector<T>>);
template <typename T>
void makeUnique(T&);
template <typename T>
std::vector<uint> getRange(std::vector<T>);
bool readLunarCraterEntry(  const std::string, 
                            lunar_crater&,
                            const char, 
                            const float = 1.2, 
                            const float = 0.9,
                            const float = 60.0,
                            const float = 1000.0);
template <typename T>
void runCraterReader(const std::string,
                     std::vector<T>&,
                     const char =',',
                     const float =1.2,
                     const float =0.9,
                     const uint =200);
template <typename T>
std::tuple<float, char> stringify_lat(const T crater);
template <typename T>
std::tuple<float, char> stringify_lon(const T crater);
template <typename T>
std::string stringify_latlon(const T crater);
std::string stringify_latlon(const float, const float);
std::tuple<float, char> stringify_lat(const float);
std::tuple<float, char> stringify_lon(const float);

/**** Template definitions ****/
template <typename T>
void runCraterReader(const std::string fname,
                     std::vector<T>& craters,
                     const char sep,
                     const float max_ell,
                     const float min_arc,
                     const uint max_n_craters) {
    T db_entry;
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
                craters.push_back(db_entry);
            }
            if(count > max_n_craters) {
                break;
            }
        }
        file.close();
    }
    else
        std::cout<<"Could not open the file: "<<fname<<std::endl;

    std::cout<<craters.size() 
             <<"craters in database."
             <<std::endl<<std::endl;
    
    // uint count = 0;
    // uint n_print = 6;
    // for (const auto& crater : entries) {
    //     printLunarCratersInfo(crater);
    //     std::cout<<std::endl;
    //     if(++count >= n_print)
    //         break;
    // }
}