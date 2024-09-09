#include "gtest/gtest.h"
#include <string>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <opencv2/core/core.hpp> 
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/viz/types.hpp>
#include <random>

#include "io.h"

std::string stringify_lat(const double lat) {
  if(lat > 90. || lat < -90) {
    std::cerr << "Latitude falls outside [-90°,90°]: " << lat << std::endl;
    return "";
  }
  std::string n_or_s;
  double abs_lat;
  n_or_s = (lat > 0)?"N":"S";
  abs_lat = abs(lat);
  std::stringstream outStream;
  outStream << std::fixed << std::setw(6) << std::fixed << std::setprecision(2) << abs_lat << "°" << n_or_s;
  return outStream.str();
}

std::string stringify_lon(const double lon) {
  double re_lon = lon;
  if(re_lon > 180)
    re_lon -= 360;
  if(re_lon > 180. || re_lon < -180) {
    std::cerr << "Longitude falls outside [-180°,180°]: " << lon << std::endl;
    return "";
  }
  double abs_lon = abs(re_lon);
  std::string e_or_w;
  e_or_w = (re_lon < 0)?"W":"E";
  std::stringstream outStream;
  outStream << std::fixed << std::setw(6) << std::fixed << std::setprecision(2) << abs_lon << "°" << e_or_w;
  return outStream.str();
}

std::string stringify_latlon(const double lat, const double lon) {
  std::string str_lat, str_lon;
  str_lat = stringify_lat(lat);
  str_lon = stringify_lon(lon); 
  std::cout << "Lat: " << lat << " | Lon: " << lon << std::endl; 
  std::cout << "Lat: " << str_lat << " | Lon: " << str_lon << std::endl; 
  return (str_lat + " " + str_lon);
}

// std::string stringify_lat(const double lat) {
//   if(lat > 90. || lat < -90) {
//     std::cerr << "Latitude falls outside [-90°,90°]: " << lat << std::endl;
//     return "";
//   }
//   std::string n_or_s;
//   double abs_lat;
//   n_or_s = (lat > 0)?"N":"S";
//   abs_lat = abs(lat);
//   std::stringstream outStream;
//   outStream << std::fixed << std::setw(6) << std::fixed << std::setprecision(2) << abs_lat << "°" << n_or_s;
//   return outStream.str();
// }

// std::string stringify_lon(const double lon) {
//   double re_lon = lon;
//   if(re_lon > 180)
//     re_lon -= 360;
//   if(re_lon > 180. || re_lon < -180) {
//     std::cerr << "Longitude falls outside [-180°,180°]: " << lon << std::endl;
//     return "";
//   }
//   double abs_lon = abs(re_lon);
//   std::string e_or_w;
//   e_or_w = (re_lon < 0)?"W":"E";
//   std::stringstream outStream;
//   outStream << std::fixed << std::setw(6) << std::fixed << std::setprecision(2) << abs_lon << "°" << e_or_w;
//   return outStream.str();
// }

// std::string stringify_latlon(const double lat, const double lon) {
//   std::string str_lat, str_lon;
//   str_lat = stringify_lat(lat);
//   str_lon = stringify_lon(lon);  
//   return (str_lat + " " + str_lon);
// }

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

class IOTest : public testing::Test {
  protected:
    IOTest() {
      // std::array<double, 5> arr = {10.0, 7.0, 300.0, 50.0, 0.0};
      // Conic conic_arr(arr);
      // c_arr.SetGeometricParameters(arr);
    }
  public:
      // // ~QueueTest() override = default;
  // Queue<int> q0_;
  // Queue<int> q1_;
  // Queue<int> q2_;
};

bool readLunarCraterEntry(lunar_crater& crater,
                          const std::string entry, 
                          const char sep,
                          const double max_ell,
                          const double min_arc,
                          const double min_diam,
                          const double max_diam) {
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
  // std::cerr << "Crater: " << crater << std::endl;
  return crater.ell < max_ell 
      && crater.arc > min_arc
      && crater.diam < max_diam
      && crater.diam > min_diam;
}

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
}


TEST(IOTest, ReadDatabaseEntry) {

  std::string string_entry = "CRATER_ID,0,0,0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0.95,0.99,500";
  char sep=',';
  double max_ell=1.1;
  double min_arc=0.9;
  double min_diam=20.;
  double max_diam=200;
  lunar_crater crater;
  bool success = io::readLunarCraterEntry(crater,
                                          string_entry, 
                                          sep,
                                          max_ell,
                                          min_arc,
                                          min_diam,
                                          max_diam);
  ASSERT_TRUE(success);                        
}

TEST(IOTest, ReadDatabase) {
  std::string fname;
  // Crater Reading
  // fname = "/home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv";
  // fname = "/home/dqr0509/data/craters/lunar_craters.csv";
  fname = "/home/ndvr/data/craters/lunar_crater_database_red2000.csv";

  char sep=',';
  double max_ell=1.2;
  double min_arc=0.9;
  double min_diam=50.;
  double max_diam=60.;
  std::vector<lunar_crater> craters;
  lunar_crater crater;                                    
  runCraterReader(craters, fname, sep, max_ell, min_arc, min_diam, max_diam, 200);

  uint c_idx = 0;
  for(auto& crater : craters) {
    // std::cout << "Lat: " << crater.lat << std::endl;
    std::cout << "Crater " << c_idx++ << " - " << crater << std::endl;
  }
  
}

TEST(IOTest, StringifyLat) {
  double lat = 80.1;
  std::string lat_north = stringify_lat(lat);
  std::string lat_cmp_n = " 80.10°N";
  int str_cmp_n = lat_north.compare(lat_cmp_n);
  EXPECT_EQ(str_cmp_n, 0);
  EXPECT_EQ(lat_north, lat_cmp_n);

  std::string lat_south = stringify_lat(-lat);
  std::string lat_cmp_s = " 80.10°S";
  int str_cmp_s = lat_south.compare(lat_cmp_s);
  EXPECT_EQ(str_cmp_s, 0);
  EXPECT_EQ(lat_south, lat_cmp_s);
  
  double lat_oob = 91;
  std::string str_oob = stringify_lat(lat_oob);  
  ASSERT_TRUE(str_oob.empty());
}

TEST(IOTest, StringifyLon) {
  double lon = 110;
  std::string lon_east = stringify_lon(lon);
  std::string lon_cmp_e = "110.00°E";
  int str_cmp_e = lon_east.compare(lon_cmp_e);
  EXPECT_EQ(str_cmp_e, 0);
  EXPECT_EQ(lon_east, lon_cmp_e);

  std::string lon_west= stringify_lon(-lon);
  std::string lon_cmp_w = "110.00°W";
  int str_cmp_w = lon_west.compare(lon_cmp_w);
  EXPECT_EQ(str_cmp_w, 0);
  EXPECT_EQ(lon_west, lon_cmp_w);
  
  double lon_oob = -181;
  std::string str_oob = stringify_lon(lon_oob);  
  EXPECT_TRUE(str_oob.empty());

  double lon_dec = 123.45;
  std::string lon_d = stringify_lon(lon_dec);
  std::string lon_d_cmp = "123.45°E";
  int str_cmp_d = lon_d.compare(lon_d_cmp);
  EXPECT_EQ(str_cmp_d, 0);
  EXPECT_EQ(lon_d, lon_d_cmp);
}

TEST(IOTest, StringifyLatLon) {
  double lat = 45.15;
  double lon = 30.95;
  std::string s_latlon = stringify_latlon(lat, lon);
  std::string s_latlon_cmp = " 45.15°N  30.95°E";
  int str_cmp = s_latlon.compare(s_latlon_cmp);
  EXPECT_EQ(str_cmp, 0);
  EXPECT_EQ(s_latlon, s_latlon_cmp);
}

