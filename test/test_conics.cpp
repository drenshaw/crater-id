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
      smajor = 100., sminor = 70., xcen = 300., ycen = 50., angle = 0.;
      arr = {smajor, sminor, xcen, ycen, angle};
      conic_a = new Conic(smajor, sminor, xcen, ycen, angle);
      conic_arr = new Conic(arr);
      // c_arr.setGeometricParameters(arr);
    }
  public:
    Conic* conic_a;
    Conic* conic_b;
    Conic* conic_arr;
    double smajor, sminor, xcen, ycen, angle;
    std::array<double, GEOMETRIC_PARAM> arr;
    // ~QueueTest() override = default;
};

TEST_F(ConicTest, CheckId) {
  // TODO: Can we ensure that this test runs first of all the conic tests?
  // TODO: Maybe the test should only ensure ID increments from one to another
  EXPECT_EQ(conic_a  ->getID(), 0);
  EXPECT_EQ(conic_arr->getID(), 1);
}

TEST_F(ConicTest, ConicInit) {
  ASSERT_DOUBLE_EQ(conic_arr->getSemiMajorAxis(), smajor);
  ASSERT_DOUBLE_EQ(conic_arr->getSemiMinorAxis(), sminor);
  ASSERT_DOUBLE_EQ(conic_arr->getCenterX(), xcen);
  ASSERT_DOUBLE_EQ(conic_arr->getCenterY(), ycen);
  ASSERT_DOUBLE_EQ(conic_arr->getAngle(), angle);
}

TEST_F(ConicTest, ConicEqual) {
  std::vector<double> vec = {smajor, sminor, xcen, ycen, angle};
  std::array<double, GEOMETRIC_PARAM> tmp_arr;
  io::copy_vec2array(vec, tmp_arr);

  Conic conic_vec(vec);
  // Conic conic_arr(tmp_arr);
  // Conic conic_var(smajor, sminor, xcen, ycen, angle);  
  EXPECT_EQ(*conic_a,   conic_arr);
  EXPECT_EQ(conic_vec, conic_arr);
}

TEST_F(ConicTest, ConicNotEqual) {
  Conic conic_major(smajor+1, sminor, xcen, ycen, angle); 
  Conic conic_minor(smajor, sminor+1, xcen, ycen, angle); 
  Conic conic_xcenter(smajor, sminor, xcen+1, ycen, angle); 
  Conic conic_ycenter(smajor, sminor, xcen, ycen+1, angle); 
  Conic conic_angle(smajor, sminor, xcen, ycen, angle+1); 

  // Messing with malloc'd objects makes things interesting
  ASSERT_NE(*conic_a, conic_major);
  ASSERT_NE(*conic_a, conic_minor);
  ASSERT_NE(*conic_a, conic_xcenter);
  ASSERT_NE(*conic_a, conic_ycenter);
  ASSERT_NE(*conic_a, conic_angle);
  ASSERT_NE(conic_major,   *conic_a);
  ASSERT_NE(conic_minor,   *conic_a);
  ASSERT_NE(conic_xcenter, *conic_a);
  ASSERT_NE(conic_ycenter, *conic_a);
  ASSERT_NE(conic_angle,   *conic_a);
}

TEST_F(ConicTest, ConicSetIndividualParameters) {
  Conic conic_empty(smajor, sminor, xcen, ycen, angle);
  // Conic conic_empty();

  // Changing values
  double d_smajor = 200., d_sminor = 60., d_xcen = 200., d_ycen = 70., d_angle = 10.;
  conic_empty.setSemimajorAxis(d_smajor);
  conic_empty.setSemiminorAxis(d_sminor);
  conic_empty.setCenterX(d_xcen);
  conic_empty.setCenterY(d_ycen);
  conic_empty.setAngle(d_angle);

  EXPECT_DOUBLE_EQ(conic_empty.getSemiMajorAxis(), d_smajor);
  EXPECT_DOUBLE_EQ(conic_empty.getSemiMinorAxis(), d_sminor);
  EXPECT_DOUBLE_EQ(conic_empty.getCenterX(), d_xcen);
  EXPECT_DOUBLE_EQ(conic_empty.getCenterY(), d_ycen);
  EXPECT_DOUBLE_EQ(conic_empty.getAngle(), d_angle);
}

TEST_F(ConicTest, ConicSetImplicit) {
  std::array<double, IMPLICIT_PARAM> impl;
  impl = conic_a->getImplicit();

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

  Eigen::Matrix3d locus = conic_a->getLocus();

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

TEST_F(ConicTest, LocusEnvelopeConversion) {
  std::array<double, IMPLICIT_PARAM> impl = {1, 0, 3, 4, 5, 7};
  double impl_norm = vectorNorm(impl);
  ASSERT_NE(impl_norm, 0);
  double last_impl = 1/impl.at(5);
  
  Conic conic_env_conversion;
  conic_env_conversion.setImplicitParameters(impl);
  normalizeImplicitParameters(impl);

  EXPECT_DOUBLE_EQ(impl.at(0), 1.*last_impl);
  EXPECT_DOUBLE_EQ(impl.at(1), 0.*last_impl);
  EXPECT_DOUBLE_EQ(impl.at(2), 3.*last_impl);
  EXPECT_DOUBLE_EQ(impl.at(3), 4.*last_impl);
  EXPECT_DOUBLE_EQ(impl.at(4), 5.*last_impl);
  EXPECT_DOUBLE_EQ(impl.at(5), 1.0);
}

TEST_F(ConicTest, ConvertEigenVectorToVector) {
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

TEST(VisualTest, EllipsePlottedCorrectly) {
  // std::vector<Conic> conics;
  // double a = 50, b = 30, r = 300;
  // int im_row = 1080, im_col = 1920;
  // for(int i = 0; i < 15; i++) {
  //   double theta = deg2rad(i*24);
  //   double xc = double(im_col)/2 + r*std::cos(theta);
  //   double yc = double(im_row)/2 + r*std::sin(theta);
  //   Conic conic(a, b, xc, yc, theta);
  //   conics.push_back(conic);
  // }

  // cv::Mat image(im_row, im_col, CV_8UC3, 
  //               cv::Scalar(50, 50, 50));
  // viz::drawEllipses(image, conics, viz::CV_colors);
  // // cv::Mat outImg;
  // // double scaling = 0.4;
  // // cv::resize(image, outImg, cv::Size(), scaling, scaling);
  // cv::imshow("Conics arrayed around center of camera", image); 
  // cv::waitKey(0); 
}


