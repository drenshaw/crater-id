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
  // TODO: Can we ensure that this test runs first of all the conic tests?
  // TODO: Maybe the test should only ensure ID increments from one to another
  Conic conic1(100.0, 70.0, 300.0, 50.0, 0.0);
  EXPECT_EQ(conic1.GetID(), 0);
  Conic conic2(100.0, 70.0, 300.0, 50.0, 0.0);
  EXPECT_EQ(conic2.GetID(), 1);
}

TEST(ConicTest, conic_init) {
  std::vector<double> vec = {100.0, 70.0, 300.0, 50.0, 0.0};
  std::array<double, GEOMETRIC_PARAM> arr;
  copy_vec2array(vec, arr);

  // TODO: make a constructor to deal with vector
  // Conic conic_vec(vec);
  Conic conic_arr(arr);
  Conic conic_var(100.0, 70.0, 300.0, 50.0, 0.0);
  Conic conic_empty;
  double semimajor_init_val = 1.0;
  
  ASSERT_EQ(conic_arr, conic_var);
  ASSERT_EQ(conic_empty.GetSemiMajorAxis(), semimajor_init_val);
}




int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}