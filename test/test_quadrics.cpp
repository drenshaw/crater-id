#include "gtest/gtest.h"
// #include <eigen3/Eigen/Dense>
// #include <random>

#include "quadrics.h"
#include "math_utils.h"

// #include <iomanip>

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
  }

public:
  // // ~QueueTest() override = default;
  // Queue<int> q0_;
  // Queue<int> q1_;
  // Queue<int> q2_;
};

TEST(QuadricTest, MakingQuadric) {
  double lat = 0, lon = 0, radius = 200;
  std::string id = "testLatLon";
  Quadric quad(lat, lon, radius, id);
  std::cerr << quad << std::endl;
  // std::cerr << quad.GetLocus() << std::endl;
  Eigen::Matrix3d t_e2m;
  // For clarity: we use a passive transformation, not an active rotation
  //   local x col vector is component in the Y-direction
  //   local y col vector is component in the Z-direction
  //   local z col vector is component in the X-direction
  t_e2m <<  0, 0, 1, 
            1, 0, 0, 
            0, 1, 0;
  Eigen::Vector3d X = Eigen::Vector3d::UnitX();
  Eigen::Vector3d Y = Eigen::Vector3d::UnitY();
  Eigen::Vector3d Z = Eigen::Vector3d::UnitZ();
  ASSERT_TRUE(quad.GetQuadricTransformationMatrix()==t_e2m);
  std::cout << "X vector pointing here : " << (t_e2m*X).transpose() << std::endl;
  std::cout << "Y vector pointing here : " << (t_e2m*Y).transpose() << std::endl;
  std::cout << "Z vector pointing here : " << (t_e2m*Z).transpose() << std::endl;
  
  Quadric north_pole(89, 0, radius, "northPole");
  std::cerr << north_pole << std::endl;
  // std::cerr << north_pole.GetLocus() << std::endl;
  std::cerr << north_pole.GetQuadricTransformationMatrix() << std::endl;
  // std::cerr << "Center to crater radius: " << calculateCraterRimFromRadius(500) << std::endl;
}

TEST(QuadricTest, InvalidCraterRadius) {
  double lat = 0, lon = 0, invalid_radius = 17000;
  std::string id_excess = "excessRadius";
  EXPECT_THROW(Quadric quad(lat, lon, invalid_radius, id_excess), std::runtime_error);
  EXPECT_THROW(Quadric quad(lat, lon, 0, "zeroRadius"), std::runtime_error);
}

TEST(QuadricTest, InvalidPositionWithDependentNormal) {
  double radius = 100;
  Eigen::Vector3d position;
  position << 0, 0, 0;
  std::string id = "dependentNormal";
  EXPECT_THROW(Quadric quad(position, radius, id), std::runtime_error);
}

TEST(QuadricTest, InvalidSurfaceNormal) {
  double lat = 0, lon = 0, radius = 17000;
  Eigen::Vector3d position = 1e5*latlon2bearing(lat, lon);
  Eigen::Vector3d surface_normal;
  surface_normal << 0, 0, 0;
  std::string id = "zeroNormal";
  EXPECT_THROW(Quadric quad(position, radius, surface_normal, id), std::runtime_error);
}
