#include "gtest/gtest.h"
#include <vector>
#include <tuple>
#include <iostream>
#include <random>

#include "structs.h"
#include "io.h"
#include "combinatorics.h"

class ComboTest : public testing::Test {
  protected:
    ComboTest() {
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

TEST(ComboTest, CheckId) {
  std::string fname;
  // Crater Reading
  // cout<<"Enter crater file name: ";
  // cin>>fname;
  fname = "/home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv";
  // fname = "/home/dqr0509/data/craters/lunar_craters.csv";

  std::vector<lunar_crater> craters;
  // runCraterReader(fname, craters, ',', 1.2, 0.9, 200);
  // uint c_idx = 0;
  // for(auto& crater : craters) {
  //   std::cout << "Crater " << c_idx++ << " " << crater << std::endl;
  // }
  
  std::vector<std::tuple<uint, uint>> valids;
  const float max_angle = 10.;
  specialCombination(craters, valids, max_angle);
  
  std::vector<std::tuple<uint, uint, uint>> triads;
  formTriads(valids, triads);

  std::cout << triads.size() << " valid triads found." << std::endl;
  // uint idx = 0;
  // for(const auto& [i, j, k]: triads) {
  //   std::cout << "IDX: " << idx++ << " | " << i << ", " << j << ", " << k << std::endl;
  // }

  const int max_iter = 10;
  print_triads(triads, craters, max_iter);
  
}

