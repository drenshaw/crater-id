#include "gtest/gtest.h"
#include <array>
#include <iostream>
#include <math.h>
#include <eigen3/Eigen/Dense>

#include "vector_math.h"

TEST(MathTest, deg2rad) {
  ASSERT_EQ(deg2rad(0), 0.0) << "Basic degrees to radians failed; verify.";
  ASSERT_LE(abs(deg2rad( 180) - M_PI), 0.000005) << "Basic +180 degrees to radians failed; verify.";
  ASSERT_LE(abs(deg2rad(-180) + M_PI), 0.000005) << "Basic -180 degrees to radians failed; verify.";
  ASSERT_DOUBLE_EQ(deg2rad( 180), M_PI);
  ASSERT_DOUBLE_EQ(deg2rad(-180), -M_PI);
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
  ASSERT_DOUBLE_EQ(vec(0), 3.);
  ASSERT_DOUBLE_EQ(vec(1), 6.);
}

TEST(MathTest, EigenDivScalar) {
  Eigen::Vector2d vec;
  vec << 6,9;
  vec /= 3;
  ASSERT_DOUBLE_EQ(vec(0), 2.);
  ASSERT_DOUBLE_EQ(vec(1), 3.);
}

TEST(MathTest, LatLon2Bearing) {
  // TODO: add more non-90 degree values for lat, lon
  double lat1, lon1 = 0;
  double epsilon = 1e-10;
  Eigen::Vector3d bearing1 = latlon2bearing(lat1, lon1);
  EXPECT_NEAR(bearing1(0), 1, epsilon);
  EXPECT_NEAR(bearing1(1), 0, epsilon);
  EXPECT_NEAR(bearing1(2), 0, epsilon);

  double lat2 = 90, lon2 = 0;
  Eigen::Vector3d bearing2 = latlon2bearing(lat2, lon2);
  EXPECT_NEAR(bearing2(0), 0, epsilon);
  EXPECT_NEAR(bearing2(1), 0, epsilon);
  EXPECT_NEAR(bearing2(2), 1, epsilon);

  double lat3 = 0, lon3 = 90;
  Eigen::Vector3d bearing3 = latlon2bearing(lat3, lon3);
  EXPECT_NEAR(bearing3(0), 0, epsilon);
  EXPECT_NEAR(bearing3(1), 1, epsilon);
  EXPECT_NEAR(bearing3(2), 0, epsilon);

  double lat4 = 0, lon4 = 180;
  Eigen::Vector3d bearing4 = latlon2bearing(lat4, lon4);
  EXPECT_NEAR(bearing4(0),-1, epsilon);
  EXPECT_NEAR(bearing4(1), 0, epsilon);
  EXPECT_NEAR(bearing4(2), 0, epsilon);

  double lat5 = -90, lon5 = 0;
  Eigen::Vector3d bearing5 = latlon2bearing(lat5, lon5);
  EXPECT_NEAR(bearing5(0), 0, epsilon);
  EXPECT_NEAR(bearing5(1), 0, epsilon);
  EXPECT_NEAR(bearing5(2),-1, epsilon);
}

TEST(MathTest, VectorNorm) {
  std::vector<double> vec1 = {3, 4};
  double norm1 = vectorNorm(vec1);
  ASSERT_DOUBLE_EQ(norm1, 5);

  std::vector<double> vec2 = {1,2,3,4,5};
  double norm2 = vectorNorm(vec2);
  double sum2 = 0;
  for(auto& element : vec2) {
      sum2 += element * element;
  }
  double norm2_check = sqrt(sum2);
  ASSERT_DOUBLE_EQ(norm2, norm2_check);

  // Checking overloaded function given iterators of vector
  double norm2_2 = vectorNorm(vec2.begin(), vec2.end());
  ASSERT_DOUBLE_EQ(norm2_2, norm2_check);
  double norm2_3 = sqrt(std::inner_product(vec2.begin(), vec2.end(), vec2.begin(), 0));
  ASSERT_DOUBLE_EQ(norm2_3, norm2_check);

  // Can we do this for std::array? Yes!
  std::array<double, 5> arr = {1,2,3,4,5};
  double norm3 = vectorNorm(arr);
  ASSERT_DOUBLE_EQ(norm3, norm2_check);

  // Can we do this for Eigen::VectorXd?
  Eigen::Vector2d evec;
  evec << 3,4;
  double norm4 = evec.norm();
  // long double norm4 = vectorNorm(evec);
  ASSERT_DOUBLE_EQ(norm4, 5);
}

TEST(MathTest, VectorCopyToStdArray) {
  const size_t n_elem = 10;
  const double val = 1.5;
  std::vector<double> vec(n_elem, val);
  std::array<double, n_elem> arr;
  copy_vec2array(vec, arr);
  for(auto& element : arr) {
    ASSERT_DOUBLE_EQ(element, val);
  }
  Eigen::Vector3d x, y;
  x << 1,2,3;
  y << 1,0,0;
  double z = x.dot(y);
  ASSERT_DOUBLE_EQ(z, 0);
}


