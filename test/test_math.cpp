#include "gtest/gtest.h"
#include <array>
#include <iostream>
#include <math.h>
#include <eigen3/Eigen/Dense>
#include <random>
#include <ctime>
// #include <unordered_set>

#include "vector_math.h"
double epsilon = 1e-15;

// template <typename T>
double random_between(const double lower_limit, const double upper_limit) {
  std::uniform_real_distribution<double> unif(lower_limit, upper_limit);
  std::default_random_engine re(time(0));
  return unif(re);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
  os << "\n[";
  for(const auto& element : vec) {
    os << "\n  " << element;
  }
  os << "\n]";
  return os;
}

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
  Eigen::Vector3d np2;
  GetNorthPoleUnitVector(np2);
  ASSERT_EQ(np2(0), 0.);
  ASSERT_EQ(np2(1), 0.);
  ASSERT_EQ(np2(2), 1.);
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
  ASSERT_DOUBLE_EQ(z, 1.0);
}

TEST(MathTest, ArgMin) {
  // Testing both std::vector and std::array
  int min_idx = 2;
  std::vector<double> vec = {3,5,1,2,4};
  std::array<double, 5> arr;
  copy_vec2array(vec, arr);
  int argmin_vec = arg_min(vec);
  int argmin_arr = arg_min(arr);
  ASSERT_EQ(argmin_vec, min_idx);
  ASSERT_EQ(argmin_arr, min_idx);

}

TEST(MathTest, ArgMax) {
  // Testing both std::vector and std::array
  int max_idx = 1;
  std::vector<double> vec = {3,5,1,2,4};
  std::array<double, 5> arr;
  copy_vec2array(vec, arr);
  int argmax_vec = arg_max(vec);
  int argmax_arr = arg_max(arr);
  ASSERT_EQ(argmax_vec, max_idx);
  ASSERT_EQ(argmax_arr, max_idx);
}

TEST(MathTest, Cofactor) {
  // TODO: fill this in later
}

TEST(MathTest, CofactorMatrix) {
  Eigen::Matrix4d mtx_4x4, cof_4x4, cof_4x4_check;
  mtx_4x4 << 5,  -2,  2,  7,
          1,   0,  0,  3,
          -3,  1,  5,  0,
          3,  -1, -9,  4;
  cof_4x4_check <<  -12, -56, 4, 4,
                    76, 208, 4, 4,
                    -60, -82, -2, 20,
                    -36, -58, -10, 12;
  cof_4x4 = getCofactorMatrix(mtx_4x4);
  ASSERT_TRUE(cof_4x4.isApprox(cof_4x4_check));

  Eigen::Matrix3d mtx_3x3, cof_3x3, cof_3x3_check;
  mtx_3x3 <<  1,9,3,
              2,5,4,
              3,7,8;
  cof_3x3_check <<  12, -4, -1,
                    -51, -1, 20,
                    21, 2, -13;
  cof_3x3 = getCofactorMatrix(mtx_3x3);
  ASSERT_TRUE(cof_3x3.isApprox(cof_3x3_check));
}

TEST(MathTest, AdjugateMatrix) {
  Eigen::Matrix4d mtx_4x4, adj_4x4, adj_4x4_check;
  mtx_4x4 <<  5, 1, -3, 3,
              1, 0, 1, -1,
              -3, 1, 5, -9,
              3, -1, -9, 4;
  adj_4x4_check <<  9, 88, -5, 4,
                    88, -224, 40, -32,
                    -5, 40, -15, -20,
                    4, -32, -20, -16;
  adj_4x4 = getAdjugateMatrix(mtx_4x4);
  ASSERT_TRUE(adj_4x4.isApprox(adj_4x4_check.transpose()));

  Eigen::Matrix3d mtx_3x3, adj_3x3, adj_3x3_check, adj_3x3_symmetric;
  mtx_3x3 <<  7, 1, -3,
              1, 4, 1,
              -3, 1, 5;
  adj_3x3_check <<  19, -8, 13,
                    -8, 26, -10,
                    13, -10, 27;
  adj_3x3 = getAdjugateMatrix(mtx_3x3);
  ASSERT_TRUE(adj_3x3.isApprox(adj_3x3_check.transpose()));

  // Check manual symmetric matrix, which has a specialized function
  adj_3x3_symmetric = get3x3SymmetricAdjugateMatrix(mtx_3x3);
  ASSERT_TRUE(adj_3x3_symmetric.isApprox(adj_3x3_check.transpose()));
}

TEST(MathTest, AngularDistance) {
  double lat1 = 0, lon1 = 0, lat2 = 0, lon2 = 90;
  Eigen::Vector3d bearing1 = latlon2bearing(lat1, lon1);
  Eigen::Vector3d bearing2 = latlon2bearing(lat2, lon2);
  double angular_distance = angularDistance(bearing1, bearing2);
  double dist_check = M_PI/2;
  ASSERT_DOUBLE_EQ(angular_distance, dist_check);
}

TEST(MathTest, AngularPseudoDistance) {
  double lat1 = 0, lon1 = 0, lat2 = 0, lon2 = 90;
  Eigen::Vector3d bearing1 = latlon2bearing(lat1, lon1);
  Eigen::Vector3d bearing2 = latlon2bearing(lat2, lon2);
  double angular_pseudodistance = angularPseudoDistance(bearing1, bearing2);
  double pseudodist_check = 0;
  // TODO: replace with using EPS; can't due to type for std::numeric_limits being undefined
  ASSERT_NEAR(angular_pseudodistance, pseudodist_check, 1e-15);
}

TEST(MathTest, LatLonDist) {
  double lat1 = 0, lon1 = 0, lat2 = 0, lon2 = 90;
  double dist = latlon_dist(lat1, lon1, lat2, lon2);
  double dist_check = M_PI/2;
  // TODO: replace with using EPS; can't due to type for std::numeric_limits being undefined
  ASSERT_NEAR(dist, dist_check, 1e-15);
}

TEST(MathTest, LLARtoECEF) {
  // TODO: Check which one is more accurate; prefer LLHtoECEF for the time being
  double radius = 6378137.0;
  double flattening = 1.0/298.257223563;

  double lat1 = 0, lon1 = 0;
  double alt1 = 0;
  Eigen::Vector3d llar1 = llarToECEF(lat1, lon1, alt1, radius, flattening);
  Eigen::Vector3d llh1  = LLHtoECEF(lat1, lon1, alt1);
  ASSERT_TRUE(llar1.isApprox(llh1));
  ASSERT_NEAR(llar1(0), radius, epsilon);
  ASSERT_NEAR(llar1(1), 0.0, epsilon);
  ASSERT_NEAR(llar1(2), 0.0, epsilon);
  // std::cerr << "LLAR (0,0): \n" << llar1 << std::endl;

  double lat2 = 37, lon2 = 21, alt2 = 500;
  Eigen::Vector3d llar2 = llarToECEF(lat2, lon2, alt2, radius, flattening);
  Eigen::Vector3d llh2  = LLHtoECEF(  lat2, lon2, alt2);
  ASSERT_TRUE(llar2.isApprox(llh2, 5e4));

  double lat3 = random_between(-90.,90.);
  double lon3 = random_between(0., 360.);
  double alt3 = random_between(0.,1000.);
  Eigen::Vector3d llar3 = llarToECEF(lat3, lon3, alt3, radius, flattening);
  Eigen::Vector3d llh3  = LLHtoECEF(  lat3, lon3, alt3);
  // std::cerr << "Lat: " << lat3 << " | " 
  //           << "Lon: " << lon3 << " | " 
  //           << "Alt: " << alt3 << std::endl;
  // std::cerr << "LLAR: \n" << llar3 << std::endl;
  // std::cerr << "LLH:  \n" << llh3  << std::endl;
  // std::cerr << "LLAR - LLH:  \n" << llar3 - llh3  << std::endl;
  // TODO: sometimes these two similar functions produce the negative of the other
  ASSERT_TRUE(llar3.isApprox(llh3, 5e4));


}

TEST(MathTest, GetRange) {
  const size_t n_elem = 4;
  const double val = 1.5;
  std::vector<double> vec(n_elem, val);
  std::array<double, n_elem> arr;

  std::vector<uint> range_vec = getRange(vec);
  ASSERT_EQ(range_vec.at(0), 0);
  ASSERT_EQ(range_vec.at(1), 1);
  ASSERT_EQ(range_vec.at(2), 2);
  ASSERT_EQ(range_vec.at(3), 3);
  double rand_val = random_between(-1.0, 1.0);
}



// template <typename T>
// void makeUniqueSet(T& vec) {
//   std::unordered_set<T> s;
//   for(const auto& element : vec) {
//     s.insert(element);
//   }
//   vec.assign(s.begin(), s.end());
//   sort(vec.begin(), vec.end());
// }

TEST(MathTest, MakeUnique) {
  std::vector<int> vec = {3,5,1,2,4,12,4,3,4,3,4,3,4};
  std::vector<int> vec_check = {1,2,3,4,5,12};
  // std::cerr << "Vector: " << vec << std::endl;
  makeUnique(vec);
  // std::cerr << "Vector: " << vec << std::endl;
  ASSERT_EQ(vec.at(0), vec_check.at(0));
  ASSERT_EQ(vec.at(1), vec_check.at(1));
  ASSERT_EQ(vec.at(2), vec_check.at(2));
  ASSERT_EQ(vec.at(3), vec_check.at(3));
  ASSERT_EQ(vec.at(4), vec_check.at(4));

}

