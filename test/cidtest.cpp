#include "gtest/gtest.h"
#include <array>
#include <iostream>

#include "conics.h"
#include "vector_math.h"
/* // CATCH2 Testing Framework
#include <catch2/catch_test_macros.hpp>

TEST_CASE( "Factorial of 0 is 1 (fail)", "[single-file]" ) {
    REQUIRE( 1 == 1 );
    std::cout << "Test1\n";
}
*/

TEST(MathTest, deg2rad) {
  ASSERT_EQ(deg2rad(0), 0.0) << "Basic degrees to radians failed; verify.";
  ASSERT_LE(abs(deg2rad( 180) - 3.141592654), 0.000005) << "Basic +180 degrees to radians failed; verify.";
  ASSERT_LE(abs(deg2rad(-180) + 3.141592654), 0.000005) << "Basic -180 degrees to radians failed; verify.";

}

TEST(MathTest, NorthPole) {
  Eigen::Vector3d np = GetNorthPoleUnitVector();
  ASSERT_EQ(np(0), 0.);
  ASSERT_EQ(np(1), 0.);
  ASSERT_EQ(np(2), 1.);
}

TEST(MathTest, EigenMultScalar) {
  Eigen::Vector2d vec;
  vec << 1,2;
  vec *= 3;
  ASSERT_EQ(vec(0), 3);
  ASSERT_EQ(vec(1), 6);
}

class ConicTest : public testing::Test {
  protected:
    ConicTest() {
      std::array<double, 5> arr = {10.0, 7.0, 300.0, 50.0, 0.0};
      // Conic conic_arr(arr);
      // c_arr.SetGeometricParameters(arr);
    }
  public:
      // // ~QueueTest() override = default;
  // Queue<int> q0_;
  // Queue<int> q1_;
  // Queue<int> q2_;
};

TEST(ConicTest, CheckId) {
  Conic conic1(100.0, 70.0, 300.0, 50.0, 0.0);
  EXPECT_EQ(conic1.GetID(), 0);
  Conic conic2(100.0, 70.0, 300.0, 50.0, 0.0);
  EXPECT_EQ(conic2.GetID(), 1);
}

TEST(ConicTest, conic_init) {
  std::array<double, 5> arr = {100.0, 70.0, 300.0, 50.0, 0.0};
  Conic conic_arr(arr);
  Conic conic_var(100.0, 70.0, 300.0, 50.0, 0.0);
  
  // ASSERT_EQ(conic_arr, conic_var);
}




int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}