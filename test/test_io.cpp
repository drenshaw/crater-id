#include "gtest/gtest.h"
#include <string>
#include <vector>
#include <array>

#include "io.h"
#include "math_utils.h"

template <typename T>
void stringifyVectorOfVectors(const std::vector<std::vector<T>> vec) {
  std::cout << "Printing vector of vectors: " << std::endl;
  for(auto& combo : vec) {
    io::stringifyVector(combo);
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



TEST(MathTest, VectorCopyToStdArray) {
  const size_t n_elem = 10;
  const double val = 1.5;
  std::vector<double> vec(n_elem, val);
  std::array<double, n_elem> arr;
  io::copy_vec2array(vec, arr);
  for(auto& element : arr) {
    ASSERT_DOUBLE_EQ(element, val);
  }
  Eigen::Vector3d x, y;
  x << 1,2,3;
  y << 1,0,0;
  double z = x.dot(y);
  ASSERT_DOUBLE_EQ(z, 1.0);
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

  // uint c_idx = 0;
  // for(auto& crater : craters) {
  //   // std::cout << "Lat: " << crater.lat << std::endl;
  //   std::cout << "Crater " << c_idx++ << " - " << crater << std::endl;
  // }
  
}

TEST(IOTest, StringifyLat) {
  double lat = 80.1;
  std::string lat_north = io::stringify_lat(lat);
  std::string lat_cmp_n = " 80.10°N";
  int str_cmp_n = lat_north.compare(lat_cmp_n);
  EXPECT_EQ(str_cmp_n, 0);
  EXPECT_EQ(lat_north, lat_cmp_n);

  std::string lat_south = io::stringify_lat(-lat);
  std::string lat_cmp_s = " 80.10°S";
  int str_cmp_s = lat_south.compare(lat_cmp_s);
  EXPECT_EQ(str_cmp_s, 0);
  EXPECT_EQ(lat_south, lat_cmp_s);
  
  double lat_oob = 91;
  std::string str_oob = io::stringify_lat(lat_oob);  
  ASSERT_TRUE(str_oob.empty());
}

TEST(IOTest, StringifyLon) {
  double lon = 110;
  std::string lon_east = io::stringify_lon(lon);
  std::string lon_cmp_e = "110.00°E";
  int str_cmp_e = lon_east.compare(lon_cmp_e);
  EXPECT_EQ(str_cmp_e, 0);
  EXPECT_EQ(lon_east, lon_cmp_e);

  std::string lon_west= io::stringify_lon(-lon);
  std::string lon_cmp_w = "110.00°W";
  int str_cmp_w = lon_west.compare(lon_cmp_w);
  EXPECT_EQ(str_cmp_w, 0);
  EXPECT_EQ(lon_west, lon_cmp_w);
  
  double lon_oob = -181;
  std::string str_oob = io::stringify_lon(lon_oob);  
  EXPECT_TRUE(str_oob.empty());

  double lon_dec = 123.45;
  std::string lon_d = io::stringify_lon(lon_dec);
  std::string lon_d_cmp = "123.45°E";
  int str_cmp_d = lon_d.compare(lon_d_cmp);
  EXPECT_EQ(str_cmp_d, 0);
  EXPECT_EQ(lon_d, lon_d_cmp);
}

TEST(IOTest, StringifyLatLon) {
  double lat = 45.15;
  double lon = 30.95;
  std::string s_latlon = io::stringify_latlon(lat, lon);
  std::string s_latlon_cmp = " 45.15°N  30.95°E";
  int str_cmp = s_latlon.compare(s_latlon_cmp);
  EXPECT_EQ(str_cmp, 0);
  EXPECT_EQ(s_latlon, s_latlon_cmp);
}

TEST(IOTest, stringifyVectorVector) {
  std::vector<double> vec_dbl = {1.1, 2.2, 3.3, 4.4};
  std::string res_dbl = io::stringifyVector(vec_dbl, "DOUBLE: ");
  std::string res_dbl_chk = "DOUBLE: 1.1, 2.2, 3.3, 4.4, ";
  int res_dbl_cmp = res_dbl.compare(res_dbl_chk);
  EXPECT_EQ(res_dbl_cmp, 0);
  EXPECT_EQ(res_dbl_chk, res_dbl);
  
  std::vector<char> vec_char = {'a', 'b', 'c', 'd', 'e'};
  std::string res_char = io::stringifyVector(vec_char, "CHAR: ");
  std::string res_char_chk = "CHAR: a, b, c, d, e, ";
  int res_char_cmp = res_char.compare(res_char_chk);
  EXPECT_EQ(res_char_cmp, 0);
  EXPECT_EQ(res_char_chk, res_char);

  std::string rep_str = "Hello";
  std::vector<std::string>  vec_str(3,rep_str);
  std::string res_str = io::stringifyVector(vec_str, "STRING: ");
  std::string res_str_chk = "STRING: Hello, Hello, Hello, ";
  int res_str_cmp = res_str.compare(res_str_chk);
  EXPECT_EQ(res_str_cmp, 0);
  EXPECT_EQ(res_str_chk, res_str);
}

TEST(IOTest, stringifyVectorArray) {
  std::array<double,4 > vec_dbl = {1.1, 2.2, 3.3, 4.4};
  std::string res_dbl = io::stringifyVector(vec_dbl, "DOUBLE: ");
  std::string res_dbl_chk = "DOUBLE: 1.1, 2.2, 3.3, 4.4, ";
  int res_dbl_cmp = res_dbl.compare(res_dbl_chk);
  EXPECT_EQ(res_dbl_cmp, 0);
  EXPECT_EQ(res_dbl_chk, res_dbl);
  
  std::array<char, 5> vec_char = {'a', 'b', 'c', 'd', 'e'};
  std::string res_char = io::stringifyVector(vec_char, "CHAR: ");
  std::string res_char_chk = "CHAR: a, b, c, d, e, ";
  int res_char_cmp = res_char.compare(res_char_chk);
  EXPECT_EQ(res_char_cmp, 0);
  EXPECT_EQ(res_char_chk, res_char);

  std::string rep_str = "Hello";
  std::array<std::string, 3>  arr_str;
  std::vector<std::string> vec_str(3,rep_str);
  io::copy_vec2array(vec_str, arr_str);
  std::string res_str = io::stringifyVector(arr_str, "STRING: ");
  std::string res_str_chk = "STRING: Hello, Hello, Hello, ";
  int res_str_cmp = res_str.compare(res_str_chk);
  EXPECT_EQ(res_str_cmp, 0);
  EXPECT_EQ(res_str_chk, res_str);
}

