#include "gtest/gtest.h"
#include <array>
#include <iostream>
#include <opencv2/core/core.hpp> 
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/viz/types.hpp>
#include <random>

#include "io.h"
#include "vector_math.h"

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

TEST(IOTest, CheckId) {
  std::string fname;
  // Crater Reading
  // cout<<"Enter crater file name: ";
  // cin>>fname;
  // fname = "/home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv";
  fname = "/home/dqr0509/data/craters/lunar_craters.csv";

  std::vector<lunar_crater> craters;
  runCraterReader(fname, craters);
  // uint c_idx = 0;
  // for(auto& crater : craters) {
  //   std::cout << "Crater " << c_idx++ << " " << crater << std::endl;
  // }
  
  std::vector<std::tuple<uint, uint>> valids;
  float max_angle = 10.;
  specialCombination(craters, valids, max_angle);
  
  std::vector<std::tuple<uint, uint, uint>> triads;
  formTriads(valids, triads);

  std::cout << triads.size() << " valid triads found." << std::endl;
  // uint idx = 0;
  // for(const auto& [i, j, k]: triads) {
  //   std::cout << "IDX: " << idx++ << " | " << i << ", " << j << ", " << k << std::endl;
  // }

  int max_iter = 10;
  print_triads(triads, craters, max_iter);
  
}

