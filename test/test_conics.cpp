#include "gtest/gtest.h"
#include <array>
// #include <eigen3/Eigen/src/Geometry/Transform.h>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <opencv2/core/core.hpp> 
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
// #include <random>

#include "conics.h"
#include "io.h"
#include "math_utils.h"
#include "visuals.h"
/* // CATCH2 Testing Framework
#include <catch2/catch_test_macros.hpp>
#include <stdexcept>

TEST_CASE( "Factorial of 0 is 1 (fail)", "[single-file]" ) {
    REQUIRE( 1 == 1 );
    std::cout << "Test1\n";
}
*/

class ConicTest : public testing::Test {
  protected:
    ConicTest() {
      // std::array<double, 5> arr = {10.0, 7.0, 300.0, 50.0, 0.0};
      // Conic conic_arr(arr);
      // c_arr.setGeometricParameters(arr);
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
  EXPECT_EQ(conic1.getID(), 0);
  Conic conic2(100.0, 70.0, 300.0, 50.0, 0.0);
  EXPECT_EQ(conic2.getID(), 1);
}

TEST(ConicTest, ConicInit) {
  double smajor = 100., sminor = 70., xcen = 300., ycen = 50., angle = 0.;
  std::array<double, GEOMETRIC_PARAM> arr = {smajor, sminor, xcen, ycen, angle};

  Conic conic_arr(arr);
  ASSERT_DOUBLE_EQ(conic_arr.getSemiMajorAxis(), smajor);
  ASSERT_DOUBLE_EQ(conic_arr.getSemiMinorAxis(), sminor);
  ASSERT_DOUBLE_EQ(conic_arr.getCenterX(), xcen);
  ASSERT_DOUBLE_EQ(conic_arr.getCenterY(), ycen);
  ASSERT_DOUBLE_EQ(conic_arr.getAngle(), angle);
}

TEST(ConicTest, ConicEqual) {
  double smajor = 100., sminor = 70., xcen = 300., ycen = 50., angle = 0.;
  std::vector<double> vec = {smajor, sminor, xcen, ycen, angle};
  std::array<double, GEOMETRIC_PARAM> arr;
  io::copy_vec2array(vec, arr);

  Conic conic_vec(vec);
  Conic conic_arr(arr);
  Conic conic_var(smajor, sminor, xcen, ycen, angle);  
  ASSERT_EQ(conic_arr, conic_var);
  ASSERT_EQ(conic_arr, conic_vec);
}

TEST(ConicTest, ConicNotEqual) {
  double smajor = 100., sminor = 70., xcen = 300., ycen = 50., angle = 0.;
  Conic conic_var(smajor, sminor, xcen, ycen, angle); 
  Conic conic_major(smajor+1, sminor, xcen, ycen, angle); 
  Conic conic_minor(smajor, sminor+1, xcen, ycen, angle); 
  Conic conic_xcenter(smajor, sminor, xcen+1, ycen, angle); 
  Conic conic_ycenter(smajor, sminor, xcen, ycen+1, angle); 
  Conic conic_angle(smajor, sminor, xcen, ycen, angle+1); 

  ASSERT_NE(conic_var, conic_major);
  ASSERT_NE(conic_var, conic_minor);
  ASSERT_NE(conic_var, conic_xcenter);
  ASSERT_NE(conic_var, conic_ycenter);
  ASSERT_NE(conic_var, conic_angle);
}

TEST(ConicTest, ConicSetIndividualParameters) {
  double smajor = 100., sminor = 70., xcen = 300., ycen = 50., angle = 0.;
  Conic conic_var(smajor, sminor, xcen, ycen, angle);

  // Changing values
  double d_smajor = 200., d_sminor = 60., d_xcen = 200., d_ycen = 70., d_angle = 10.;
  conic_var.setSemimajorAxis(d_smajor);
  conic_var.setSemiminorAxis(d_sminor);
  conic_var.setCenterX(d_xcen);
  conic_var.setCenterY(d_ycen);
  conic_var.setAngle(d_angle);

  EXPECT_DOUBLE_EQ(conic_var.getSemiMajorAxis(), d_smajor);
  EXPECT_DOUBLE_EQ(conic_var.getSemiMinorAxis(), d_sminor);
  EXPECT_DOUBLE_EQ(conic_var.getCenterX(), d_xcen);
  EXPECT_DOUBLE_EQ(conic_var.getCenterY(), d_ycen);
  EXPECT_DOUBLE_EQ(conic_var.getAngle(), d_angle);
}

TEST(ConicTest, ConicSetImplicit) {
  double smajor = 100., sminor = 70., xcen = 300., ycen = 50., angle = 0.;
  Conic conic_var(smajor, sminor, xcen, ycen, angle);
  std::array<double, IMPLICIT_PARAM> impl;
  impl = conic_var.getImplicit();

  // Comparing with MATLAB version of software to ensure alignment and correct values
  // double impl_a = 1.1750273695215061e-05;
  // double impl_b = 0.0;
  // double impl_c = 2.3980150398398082e-05;
  // double impl_d = -0.0070501642171290364;
  // double impl_e = -0.0023980150398398084;
  // double impl_f = 0.99997227161320001;
  double impl_a =  1.1750599520383693e-05;
  double impl_b =  0.0;
  double impl_c =  2.3980815347721824e-05;
  double impl_d = -0.0070503597122302166;
  double impl_e = -0.0023980815347721825;
  double impl_f =  1.0;

  EXPECT_DOUBLE_EQ(impl.at(0), impl_a);
  EXPECT_DOUBLE_EQ(impl.at(1), impl_b);
  EXPECT_DOUBLE_EQ(impl.at(2), impl_c);
  EXPECT_DOUBLE_EQ(impl.at(3), impl_d);
  EXPECT_DOUBLE_EQ(impl.at(4), impl_e);
  EXPECT_DOUBLE_EQ(impl.at(5), impl_f);

  Eigen::Matrix3d locus;
  locus = conic_var.getLocus();

  EXPECT_DOUBLE_EQ(locus(0,0), impl_a);
  EXPECT_DOUBLE_EQ(locus(0,1), impl_b/2);
  EXPECT_DOUBLE_EQ(locus(1,1), impl_c);
  EXPECT_DOUBLE_EQ(locus(0,2), impl_d/2);
  EXPECT_DOUBLE_EQ(locus(1,2), impl_e/2);
  EXPECT_DOUBLE_EQ(locus(2,2), impl_f);

  EXPECT_DOUBLE_EQ(locus(0,1), locus(1,0));
  EXPECT_DOUBLE_EQ(locus(0,2), locus(2,0));
  EXPECT_DOUBLE_EQ(locus(1,2), locus(2,1));
}

TEST(ConicTest, LocusEnvelopeConversion) {
  std::array<double, IMPLICIT_PARAM> impl = {1, 0, 3, 4, 5, 7};
  double impl_norm = vectorNorm(impl);
  ASSERT_NE(impl_norm, 0);
  double last_impl = 1/impl.at(5);
  
  Conic conic_var;
  conic_var.setImplicitParameters(impl);
  normalizeImplicitParameters(impl);

  EXPECT_DOUBLE_EQ(impl.at(0), 1.*last_impl);
  EXPECT_DOUBLE_EQ(impl.at(1), 0.*last_impl);
  EXPECT_DOUBLE_EQ(impl.at(2), 3.*last_impl);
  EXPECT_DOUBLE_EQ(impl.at(3), 4.*last_impl);
  EXPECT_DOUBLE_EQ(impl.at(4), 5.*last_impl);
  EXPECT_DOUBLE_EQ(impl.at(5), 1.0);
}

TEST(ConicTest, ConvertEigenVectorToVector) {
  std::array<double, CONIC_DIM> arr;
  std::vector<double> vec;
  vec.reserve(CONIC_DIM);
  Eigen::Vector3d evec;
  evec << 1,2,3;
  convertEigenVectorToArray(evec, arr);
  convertEigenVectorToVector(evec, vec);
  EXPECT_EQ(arr.at(0), evec(0));
  EXPECT_EQ(arr.at(1), evec(1));
  EXPECT_EQ(arr.at(2), evec(2));
  EXPECT_EQ(arr.at(2), 3);
  EXPECT_EQ(vec[0], evec(0));
  EXPECT_EQ(vec[1], evec(1));
  EXPECT_EQ(vec[2], evec(2));
  EXPECT_EQ(vec[2], 3);
}

TEST(VisualTest, SlopeInterceptInvalid) {
  double A = 0, B = 0, C = 0;
  double slope, intercept;
  Eigen::Vector3d line;
  line << A, B, C;
  // viz::getSlopeInterceptFromStandard(line, slope, intercept);
  EXPECT_THROW(viz::getSlopeInterceptFromStandard(line, slope, intercept), std::runtime_error);
  ASSERT_TRUE(std::isnan(slope));
  ASSERT_TRUE(std::isnan(intercept));
}

TEST(VisualTest, SlopeInterceptZeroSlope) {
  double A = 0, B = 1, C = 5;
  double slope, intercept;
  Eigen::Vector3d line;
  line << A, B, C;
  viz::getSlopeInterceptFromStandard(line, slope, intercept);
  ASSERT_DOUBLE_EQ(slope, 0);
  ASSERT_DOUBLE_EQ(intercept, -C);
}

TEST(VisualTest, SlopeInterceptInfSlope) {
  double A = 1, B = 0, C = 5;
  double slope, intercept;
  Eigen::Vector3d line;
  line << A, B, C;
  viz::getSlopeInterceptFromStandard(line, slope, intercept);
  ASSERT_TRUE(std::isinf(slope));
  ASSERT_DOUBLE_EQ(intercept, -C);
}


