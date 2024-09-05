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

TEST(ConicTest, CheckId) {
  // TODO: Can we ensure that this test runs first of all the conic tests?
  // TODO: Maybe the test should only ensure ID increments from one to another
  Conic conic1(100.0, 70.0, 300.0, 50.0, 0.0);
  EXPECT_EQ(conic1.GetID(), 0);
  Conic conic2(100.0, 70.0, 300.0, 50.0, 0.0);
  EXPECT_EQ(conic2.GetID(), 1);
}

TEST(ConicTest, ConicInit) {
  double smajor = 100., sminor = 70., xcen = 300., ycen = 50., angle = 0.;
  std::array<double, GEOMETRIC_PARAM> arr = {smajor, sminor, xcen, ycen, angle};

  Conic conic_arr(arr);
  ASSERT_DOUBLE_EQ(conic_arr.GetSemiMajorAxis(), smajor);
  ASSERT_DOUBLE_EQ(conic_arr.GetSemiMinorAxis(), sminor);
  ASSERT_DOUBLE_EQ(conic_arr.GetCenterX(), xcen);
  ASSERT_DOUBLE_EQ(conic_arr.GetCenterY(), ycen);
  ASSERT_DOUBLE_EQ(conic_arr.GetAngle(), angle);
}

TEST(ConicTest, ConicEqual) {
  double smajor = 100., sminor = 70., xcen = 300., ycen = 50., angle = 0.;
  std::vector<double> vec = {smajor, sminor, xcen, ycen, angle};
  std::array<double, GEOMETRIC_PARAM> arr;
  copy_vec2array(vec, arr);

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
  conic_var.SetSemimajorAxis(d_smajor);
  conic_var.SetSemiminorAxis(d_sminor);
  conic_var.SetCenterX(d_xcen);
  conic_var.SetCenterY(d_ycen);
  conic_var.SetAngle(d_angle);

  EXPECT_DOUBLE_EQ(conic_var.GetSemiMajorAxis(), d_smajor);
  EXPECT_DOUBLE_EQ(conic_var.GetSemiMinorAxis(), d_sminor);
  EXPECT_DOUBLE_EQ(conic_var.GetCenterX(), d_xcen);
  EXPECT_DOUBLE_EQ(conic_var.GetCenterY(), d_ycen);
  EXPECT_DOUBLE_EQ(conic_var.GetAngle(), d_angle);
}

TEST(ConicTest, ConicSetImplicit) {
  double smajor = 100., sminor = 70., xcen = 300., ycen = 50., angle = 0.;
  Conic conic_var(smajor, sminor, xcen, ycen, angle);
  std::array<double, IMPLICIT_PARAM> impl;
  impl = conic_var.GetImplicit();

  double impl_a = 1.1750273695215061e-05;
  double impl_b = 0.0;
  double impl_c = 2.3980150398398082e-05;
  double impl_d = -0.0070501642171290364;
  double impl_e = -0.0023980150398398084;
  double impl_f = 0.99997227161320001;

  EXPECT_DOUBLE_EQ(impl.at(0), impl_a);
  EXPECT_DOUBLE_EQ(impl.at(1), impl_b);
  EXPECT_DOUBLE_EQ(impl.at(2), impl_c);
  EXPECT_DOUBLE_EQ(impl.at(3), impl_d);
  EXPECT_DOUBLE_EQ(impl.at(4), impl_e);
  EXPECT_DOUBLE_EQ(impl.at(5), impl_f);

  Eigen::Matrix3d locus;
  locus = conic_var.GetLocus();

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
  
  Conic conic_var;
  conic_var.SetImplicitParameters(impl);
  conic_var.NormalizeImplicitParameters(impl);
  EXPECT_DOUBLE_EQ(impl.at(0), 1./impl_norm);
  EXPECT_DOUBLE_EQ(impl.at(1), 0./impl_norm);
  EXPECT_DOUBLE_EQ(impl.at(2), 3./impl_norm);
  EXPECT_DOUBLE_EQ(impl.at(3), 4./impl_norm);
  EXPECT_DOUBLE_EQ(impl.at(4), 5./impl_norm);
  EXPECT_DOUBLE_EQ(impl.at(5), 7./impl_norm);
}

TEST(ConicTest, ConicIntersection) {
  double smajor = 100., sminor = 70., xcen = 300., ycen = 50., angle = 0.;
  Conic conicA(smajor, sminor, xcen, ycen, angle);
  Conic conicB(smajor, sminor, xcen-315, ycen, angle);
  // Eigen::Vector3d g; g.fill(0);
  // Eigen::Vector3d h; h.fill(0);
  std::tuple<Eigen::Vector3d, Eigen::Vector3d> gh;
  // std::get<0>(gh) = g;
  // std::get<1>(gh) = h;
  bool success = conicA.ConicIntersectionLines(conicB, gh);
  std::cerr << "Intersection: \n" << std::get<0>(gh) << std::endl << std::get<1>(gh) << std::endl;
  ASSERT_EQ(success, true);

}

