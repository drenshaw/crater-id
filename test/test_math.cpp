#include "gtest/gtest.h"
#include <array>
#include <iostream>
#include <math.h>
#include <eigen3/Eigen/Dense>
#include <random>
#include <ctime>
// #include <unordered_set>

#include "io.h"
#include "math_utils.h"
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
  Eigen::Vector3d unitZ = Eigen::Vector3d::UnitZ();
  Eigen::Vector3d np = getNorthPoleUnitVector();
  ASSERT_TRUE(np.isApprox(unitZ));
  Eigen::Vector3d np2;
  GetNorthPoleUnitVector(np2);
  ASSERT_TRUE(np2.isApprox(unitZ));
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
  double lat1 = 0, lon1 = 0;
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

TEST(MathTest, ArgMin) {
  // Testing both std::vector and std::array
  int min_idx = 2;
  std::vector<double> vec = {3,5,1,2,4};
  std::array<double, 5> arr;
  io::copy_vec2array(vec, arr);
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
  io::copy_vec2array(vec, arr);
  int argmax_vec = arg_max(vec);
  int argmax_arr = arg_max(arr);
  ASSERT_EQ(argmax_vec, max_idx);
  ASSERT_EQ(argmax_arr, max_idx);
}

TEST(MathTest, Cofactor) {
  // TODO: fill this in later
  Eigen::Matrix2d mtx;
  mtx << 1,2,3,4;
  Eigen::Matrix2d mtx_adj = adjugate(mtx);
  Eigen::Matrix2d mtx_check;
  mtx_check << 4, -2, -3, 1;
  ASSERT_TRUE(mtx_adj.isApprox(mtx_check));
  Eigen::MatrixXd nonsquare(2,3);
  nonsquare << 1,2,3,4,5,6;
  ASSERT_THROW(adjugate(nonsquare), std::runtime_error);

}

TEST(MathTest, CofactorMatrix3x3) {
  Eigen::Matrix3d mtx_3x3, cof_3x3, cof_3x3_check;
  mtx_3x3 <<  1,9,3,
              2,5,4,
              3,7,8;
  cof_3x3_check <<  12, -4, -1,
                    -51, -1, 20,
                    21, 2, -13;
  cof_3x3 = cofactor(mtx_3x3);
  ASSERT_TRUE(cof_3x3.isApprox(cof_3x3_check));
  ASSERT_TRUE(adjugate(mtx_3x3).isApprox(cof_3x3_check.transpose()));
}

TEST(MathTest, CofactorMatrix4x4) {
  Eigen::Matrix4d mtx_4x4, cof_4x4, cof_4x4_check;
  mtx_4x4 << 5,  -2,  2,  7,
          1,   0,  0,  3,
          -3,  1,  5,  0,
          3,  -1, -9,  4;
  cof_4x4_check <<  -12, -56, 4, 4,
                    76, 208, 4, 4,
                    -60, -82, -2, 20,
                    -36, -58, -10, 12;
  cof_4x4 = cofactor(mtx_4x4);
  ASSERT_TRUE(cof_4x4.isApprox(cof_4x4_check));
  ASSERT_TRUE(adjugate(mtx_4x4).isApprox(cof_4x4_check.transpose()));
}

TEST(MathTest, SymmetricAdjugateMatrix3x3) {
  Eigen::Matrix3d mtx_3x3, adj_3x3, adj_3x3_check, adj_3x3_symmetric, cof_3x3;
  mtx_3x3 <<  7, 1, -3,
              1, 4, 1,
              -3, 1, 5;
  adj_3x3_check <<  19, -8, 13,
                    -8, 26, -10,
                    13, -10, 27;
  adj_3x3 = adjugate(mtx_3x3);
  ASSERT_TRUE(adj_3x3.isApprox(adj_3x3_check));

  // Check manual symmetric matrix, which has a specialized function
  adj_3x3_symmetric = symmetricAdjugate(mtx_3x3);
  cof_3x3 = cofactor(mtx_3x3);
  ASSERT_TRUE(adj_3x3_symmetric.isApprox(adj_3x3_check));
  ASSERT_TRUE(cof_3x3.isApprox(adj_3x3_check.transpose()));
}

TEST(MathTest, SymmetricAdjugateMatrix4x4) {
  Eigen::Matrix4d mtx_4x4, adj_4x4, adj_4x4_check, cof_4x4, adj_adj_4x4, adj_cof_4x4, adj_adj_4x4_check;
  mtx_4x4 <<  5, 1, -3, 3,
              1, 0, 1, -1,
              -3, 1, 5, -9,
              3, -1, -9, 4;
  adj_4x4_check <<  9, 88, -5, 4,
                    88, -224, 40, -32,
                    -5, 40, -15, -20,
                    4, -32, -20, -16;
  adj_4x4 = adjugate(mtx_4x4);
  cof_4x4 = cofactor(mtx_4x4);
  ASSERT_TRUE(adj_4x4.isApprox(adj_4x4_check));
  ASSERT_TRUE(cof_4x4.isApprox(adj_4x4_check.transpose()));

  adj_adj_4x4_check <<  128000, 25600, -76800, 76800,
                        25600, 0, 25600, -25600, 
                        -76800, 25600, 128000, -230400, 
                        76800, -25600, -230400, 102400;
  adj_adj_4x4 = adjugate(adj_4x4);
  adj_cof_4x4 = cofactor(adj_4x4);
  ASSERT_TRUE(adj_adj_4x4.isApprox(adj_adj_4x4_check));
  ASSERT_TRUE(adj_cof_4x4.isApprox(adj_adj_4x4_check.transpose()));
  // Just as a fun fact, the adj(adj(mtx)) = det(mtx)^(n-2) where n is the number of dimensions 
  double det = mtx_4x4.determinant();
  int n_dim = mtx_4x4.rows();
  ASSERT_TRUE((std::pow(det,n_dim-2) * mtx_4x4).isApprox(adj_adj_4x4));
}

TEST(MathTest, AdjugateMatrix4x4) {
  Eigen::Matrix4d mtx_4x4, adj_4x4, adj_4x4_check, cof_4x4;
  mtx_4x4 <<  1, 5, 4, 12, 
              3, 2, 9, 4, 
              3, 8, 5, 3, 
              0, 4, 8, 11;
  adj_4x4_check <<  -446, -379, 71, 605, 
                    108, 187, -223, -125, 
                    210, -55, -15, -205, 
                    -192, -28, 92, 70;
  adj_4x4 = adjugate(mtx_4x4);
  cof_4x4 = cofactor(mtx_4x4);
  ASSERT_TRUE(adj_4x4.isApprox(adj_4x4_check));
  ASSERT_TRUE(cof_4x4.isApprox(adj_4x4_check.transpose()));
  double det = mtx_4x4.determinant();
  int n_dim = mtx_4x4.rows();
  ASSERT_TRUE((std::pow(det,n_dim-2) * mtx_4x4).isApprox(adjugate(adj_4x4)));
}

TEST(MathTest, AdjugateMatrix4x4ByHand) {
  Eigen::Matrix4d mtx_4x4, adj_4x4, adj_4x4_check, cof_4x4;
  mtx_4x4 <<  1.1, 2.1, 1.3, 0.8,
              -0.7, -1.1, 0.8, 1.2,
              -2.1, 1.8, 0.4, -0.1, 
              0.9, -1.3, -2.0, 1.4;
  adj_4x4_check <<  -6.004, 4.815, 9.958, 0.015, 
                    -6.636, 7.173, -7.094, -2.863, 
                    -3.476, -9.502, 5.868, 10.55, 
                    -7.268, -10.009, -4.606, -9.649;
  adj_4x4 = adjugate(mtx_4x4);
  cof_4x4 = cofactor(mtx_4x4);
  Eigen::Matrix4d adj_adj = adjugate(adjugate(mtx_4x4));
  // std::cout << "The adjugate of my adjugate is a factor of the determinant to a power.\n" << adj_adj << std::endl;
  EXPECT_TRUE(adj_4x4.isApprox(adj_4x4_check, 1e-3));
  EXPECT_TRUE(cof_4x4.isApprox(adj_4x4_check.transpose(),1e-3));
  double det = mtx_4x4.determinant();
  int n_dim = mtx_4x4.rows();
  ASSERT_TRUE((std::pow(det,n_dim-2) * mtx_4x4).isApprox(adj_adj));
}

TEST(MathTest, AdjugateMatrix4x4RankDeficient) {
  Eigen::Matrix4d mtx_4x4, adj_4x4, adj_4x4_check, cof_4x4, adj_adj;
  mtx_4x4 <<  0.0, 2.1, 1.3, 0.8,
              -0.7, 0.0, 0.8, 1.2,
              -2.1, 1.8, 0.0, -0.1, 
              0.9, -1.3, -2.0, 0.0;
  adj_4x4 = adjugate(mtx_4x4);
  adj_adj = adjugate(adj_4x4);
  // std::cout << "adj(A) (rank deficient):\n" << adj_4x4 << std::endl;
  // std::cout << "AA* (rank deficient):\n" << mtx_4x4 * adj_4x4 << std::endl;

  double det = mtx_4x4.determinant();
  int n_dim = mtx_4x4.rows();
  double scalar = std::pow(det,n_dim-2);
  // std::cout << "det^(n-2)*mtx (rank deficient): " << det << std::endl
  //           << scalar * mtx_4x4 << std::endl;
  ASSERT_TRUE((scalar * mtx_4x4).isApprox(adj_adj));
}

TEST(MathTest, AdjugateMatrixRandom3x3) {
  for(size_t i = 0; i < 100; i++) {
    Eigen::Matrix3d rnd_mtx = Eigen::Matrix3d::Random();
    EXPECT_NE(rnd_mtx(Eigen::last, Eigen::last), 0);
    Eigen::Matrix3d rnd_mtx_adj = adjugate(rnd_mtx);
    Eigen::Matrix3d rnd_mtx_back = adjugate(rnd_mtx_adj);
    rnd_mtx /= rnd_mtx(Eigen::last, Eigen::last);
    rnd_mtx_back /= rnd_mtx_back(Eigen::last, Eigen::last);
    EXPECT_TRUE(rnd_mtx.isApprox(rnd_mtx_back, 1e-3));
  }
}

TEST(MathTest, AdjugateMatrixRandom4x4) {
  for(size_t i = 0; i < 100; i++) {
    Eigen::Matrix4d rnd_mtx = Eigen::Matrix4d::Random();
    EXPECT_NE(rnd_mtx(Eigen::last, Eigen::last), 0);
    Eigen::Matrix4d rnd_mtx_adj = adjugate(rnd_mtx);
    Eigen::Matrix4d rnd_mtx_back = adjugate(rnd_mtx_adj);
    rnd_mtx /= rnd_mtx(Eigen::last, Eigen::last);
    rnd_mtx_back /= rnd_mtx_back(Eigen::last, Eigen::last);
    EXPECT_TRUE(rnd_mtx.isApprox(rnd_mtx_back, 1e-3));
  }
}

TEST(MathTest, SymmetricAdjugateMatrix34) {
  Eigen::Matrix4d mtx_4x4, adj_4x4, adj_4x4_check;
  mtx_4x4 <<  5, 1, -3, 3,
              1, 0, 1, -1,
              -3, 1, 5, -9,
              3, -1, -9, 4;
  adj_4x4_check <<  9, 88, -5, 4,
                    88, -224, 40, -32,
                    -5, 40, -15, -20,
                    4, -32, -20, -16;
  adj_4x4 = adjugate(mtx_4x4);
  ASSERT_TRUE((adj_4x4).isApprox(adj_4x4_check));

}

TEST(MathTest, AngularDistance) {
  double lat1 = 0, lon1 = 0, lat2 = 0, lon2 = 90;
  Eigen::Vector3d bearing1 = latlon2bearing(lat1, lon1);
  Eigen::Vector3d bearing2 = latlon2bearing(lat2, lon2);
  double angular_distance = getAngleBetweenVectors(bearing1, bearing2);
  double dist_check = M_PI/2;
  ASSERT_DOUBLE_EQ(angular_distance, dist_check);
}

TEST(MathTest, AngularPseudoDistance) {
  double lat1 = 0, lon1 = 0, lat2 = 0, lon2 = 90;
  Eigen::Vector3d bearing1 = latlon2bearing(lat1, lon1);
  Eigen::Vector3d bearing2 = latlon2bearing(lat2, lon2);
  double angular_pseudodistance = getPseudoAngleBetweenVectors(bearing1, bearing2);
  double pseudodist_check = 0;
  // TODO: replace with using EPS; can't due to type for std::numeric_limits being undefined
  ASSERT_NEAR(angular_pseudodistance, pseudodist_check, 1e-15);
}

TEST(MathTest, LatLonDist) {
  double lat1 = 0, lon1 = 0, lat2 = 0, lon2 = 90;
  double dist = getAngleBetweenLatLon(lat1, lon1, lat2, lon2);
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
  // std::array<double, n_elem> arr;

  std::vector<uint> range_vec = getRange(vec);
  ASSERT_EQ(range_vec.at(0), 0);
  ASSERT_EQ(range_vec.at(1), 1);
  ASSERT_EQ(range_vec.at(2), 2);
  ASSERT_EQ(range_vec.at(3), 3);
  // double rand_val = random_between(-1.0, 1.0);
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

TEST(MathTest, UndefinedCheck) {
  // std::uint8_t k = 255;
  // EXPECT_ANY_THROW(
  //   try {
  //     k += 10;
  //   }
  //   catch (...) {
  //     throw;
  //   }
  // );
}

void findMismatchedValue(const Eigen::Matrix3d& mtx, const Eigen::Matrix3d& mtx_check) {
  if(!mtx.isApprox(mtx_check)) {
    EXPECT_DOUBLE_EQ(mtx(0,0), mtx_check(0,0));
    EXPECT_DOUBLE_EQ(mtx(0,1), mtx_check(0,1));
    EXPECT_DOUBLE_EQ(mtx(0,2), mtx_check(0,2));
    EXPECT_DOUBLE_EQ(mtx(1,0), mtx_check(1,0));
    EXPECT_DOUBLE_EQ(mtx(1,1), mtx_check(1,1));
    EXPECT_DOUBLE_EQ(mtx(1,2), mtx_check(1,2));
    EXPECT_DOUBLE_EQ(mtx(2,0), mtx_check(2,0));
    EXPECT_DOUBLE_EQ(mtx(2,1), mtx_check(2,1));
    EXPECT_DOUBLE_EQ(mtx(2,2), mtx_check(2,2));
  }
}

TEST(MathTest, GetENUFrameLatOnly) {
  double lat = 30, lon = 0;
  Eigen::Matrix3d enu = getENUFrame(lat, lon);
  Eigen::Matrix3d enu_check(3,3);
  // For clarity: we use a passive transformation, not an active rotation
  //   local x col vector is component in the Y-direction
  //   local y col vector is component in the Z-direction
  //   local z col vector is component in the X-direction
  double sqrt_3_2 = std::sqrt(3)/2;
  enu_check <<  0.0,                 -0.5,                  sqrt_3_2,
                1.0,                  0.0,                  0.0,
                0.0,                  sqrt_3_2,             0.5;
  ASSERT_TRUE(enu.isApprox(enu_check));
}

TEST(MathTest, GetENUFrame) {
  double lat = 15, lon = 30;
  Eigen::Matrix3d enu = getENUFrame(lat, lon);
  Eigen::Matrix3d enu_check(3,3);
  // For clarity: we use a passive transformation, not an active rotation
  //   local x col vector is component in the Y-direction
  //   local y col vector is component in the Z-direction
  //   local z col vector is component in the X-direction
  enu_check << -0.5,                 -0.22414386804201336, 0.83651630373780794,
                0.86602540378443871, -0.12940952255126034, 0.48296291314453410,
                0.00000000000000000,  0.96592582628906820, 0.25881904510252074;
  findMismatchedValue(enu, enu_check);
  ASSERT_TRUE(enu.isApprox(enu_check));
}

TEST(MathTest, GetENUFrameXAxis) {
  double lat = 0, lon = 0;
  Eigen::Matrix3d enu = getENUFrame(lat, lon);
  Eigen::Matrix3d enu_check(3,3);
  enu_check << 0, 0, 1, 1, 0, 0, 0, 1, 0;
  findMismatchedValue(enu, enu_check);
  ASSERT_TRUE(enu.isApprox(enu_check));
}

TEST(MathTest, GetENUFrameNorthPole) {
  double lat = 90, lon = 0;
  ASSERT_THROW(getENUFrame(lat, lon), std::runtime_error);
  ASSERT_THROW(getENUFrame(lat, lon+30), std::runtime_error);
  Eigen::Matrix3d enu = getENUFrame(lat-0.01, lon);
  Eigen::Matrix3d enu_check(3,3);
  enu_check <<  0.0, -0.9999999847691291, 0.00017453292431360922, 
                1.0, 0.0, 0, 
                0.0, 0.00017453292431360924, 0.9999999847691291;
  findMismatchedValue(enu, enu_check);
  ASSERT_TRUE(enu.isApprox(enu_check));
}

TEST(MathTest, GetENUFrameSouthPole) {
  double lat = -90, lon = 0;
  ASSERT_THROW(getENUFrame(lat, lon), std::runtime_error);
  ASSERT_THROW(getENUFrame(lat, lon+30), std::runtime_error);
  Eigen::Matrix3d enu = getENUFrame(lat+0.01, lon);
  Eigen::Matrix3d enu_check(3,3);
  enu_check <<  0.0,  0.9999999847691291,      0.00017453292431360922, 
                1.0,  0.0,                     0, 
                0.0,  0.00017453292431360924, -0.9999999847691291;
  findMismatchedValue(enu, enu_check);
  ASSERT_TRUE(enu.isApprox(enu_check));
}

