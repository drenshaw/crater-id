#include "gtest/gtest.h"
// #include <eigen3/Eigen/Dense>
// #include <random>

#include "quadrics.h"
/* // CATCH2 Testing Framework
#include <catch2/catch_test_macros.hpp>
#include <stdexcept>

TEST_CASE( "Factorial of 0 is 1 (fail)", "[single-file]" ) {
    REQUIRE( 1 == 1 );
    std::cout << "Test1\n";
}
*/

class QuadricTest : public testing::Test {
  protected:
    QuadricTest() {
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

TEST(QuadricTest, MakingQuadric) {
  double lat = 0, lon = 0, radius = 100;
  std::string id = "testLatLon";
  Quadric quad(lat, lon, radius, id);
  std::cerr << quad << std::endl;
}
