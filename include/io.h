#ifndef IO_H
#define IO_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <tuple>
#include "structs.h"

#define NDIM 2
#define R_MOON 1737.4

std::ostream& operator<<(std::ostream& os, const lunar_crater& crater);

namespace io {

template <typename T>
void printVector(const std::vector<T>, const std::string="");
template <typename T>
void printVectorOfVectors(const std::vector<std::vector<T>>);
template <typename T>
void makeUnique(T&);
// template <typename T>
// std::vector<uint> getRange(std::vector<T>); //TODO: this exists in vector_math.h as well
bool readLunarCraterEntry(lunar_crater& crater,
                          const std::string entry, 
                          const char sep=',',
                          const double max_ell=1.2,
                          const double min_arc=0.9,
                          const double min_diam=50.,
                          const double max_diam=200.);
std::string stringify_lat(const double);
std::string stringify_lon(const double);
std::string stringify_latlon(const double, const double);
std::string stringify_lat(const lunar_crater crater);
std::string stringify_lon(const lunar_crater crater);
std::string stringify_latlon(const lunar_crater crater);

/**** Template definitions ****/
// template <typename T>
// std::string stringify_lat(const T crater) {
//   return stringify_lat(crater.lat);
// }

// template <typename T>
// std::string stringify_lon(const T crater) {
//   return stringify_lon(crater.lon);
// }

// template <typename T>
// std::string stringify_latlon(const T crater) {
//   return stringify_latlon(crater.lat, crater.lon);
// }

template <typename T>
void runCraterReader( std::vector<T>& craters,
                      const std::string fname,
                      const char sep=',',
                      const double max_ell=1.2,
                      const double min_arc=0.9,
                      const double min_diam=50.,
                      const double max_diam=200.,
                      const uint max_n_craters=200) {
  T crater;
  std::string line;

  std::ifstream file (fname, std::ios::in);
  if(file.is_open()) {
    uint count = 0;

    while(getline(file, line)) {
      std::stringstream str(line);
      if(count==0) {
        count++;
        continue;
      }
      if(readLunarCraterEntry(crater, line, sep, max_ell, min_arc, min_diam, max_diam)) {
        count++;
        craters.push_back(crater);
      }
      if(count > max_n_craters) {
        break;
      }
    }
    file.close();
  }
  else
    std::cerr<<"Could not open the file: "<<fname<<std::endl;

  std::cout << craters.size() << " craters in database.\n";
  
  // uint count = 0;
  // uint n_print = 6;
  // for (const auto& crater : entries) {
  //   printLunarCratersInfo(crater);
  //   std::cout<<std::endl;
  //   if(++count >= n_print)
  //     break;
  // }
}

template <typename T>
void printVector(const std::vector<T> vec, const std::string prepend) {
  std::cout << prepend;
  for(auto& idx : vec) {
    std::cout << idx << ", ";
  }
  std::cout << std::endl;
}

template <typename T, size_t SIZE>
void printVector(const std::array<T, SIZE> arr, const std::string prepend) {
  std::cout << prepend;
  for(auto& elem : arr) {
    std::cout << elem << ", ";
  }
  std::cout << std::endl;
}

template <typename T>
void printVectorOfVectors(const std::vector<std::vector<T>> vec) {
  std::cout << "Printing vector of vectors: " << std::endl;
  for(auto& combo : vec) {
    printVector(combo);
  }
  std::cout << std::endl;
}

} // namespace

#endif