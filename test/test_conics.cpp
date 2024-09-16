#include "gtest/gtest.h"
#include <array>
// #include <eigen3/Eigen/src/Geometry/Transform.h>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <opencv2/core/core.hpp> 
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
// #include <random>

#include "conics.h"
#include "io.h"
#include "vector_math.h"
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
  double last_impl = 1/impl.at(5);
  
  Conic conic_var;
  conic_var.SetImplicitParameters(impl);
  conic_var.NormalizeImplicitParameters(impl);

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
  convertEigenVectorToVector(evec, arr);
  convertEigenVectorToVector(evec, vec);
  EXPECT_EQ(arr.at(0), evec(0));
  EXPECT_EQ(arr.at(1), evec(1));
  EXPECT_EQ(arr.at(2), evec(2));
  EXPECT_EQ(arr.at(2), 3);
  EXPECT_EQ(vec[0], evec(0));
  EXPECT_EQ(vec[1], evec(1));
  EXPECT_EQ(vec[2], evec(2));
  EXPECT_EQ(vec[2], 3);
  // ASSERT_TRUE(success);
}

TEST(ConicTest, ConicIntersection) {
  double smajor = 70., sminor = 50., xcen = 200., ycen = 200., angle = -15.;
  Conic conicA(smajor, sminor, xcen, ycen, angle+50);
  Conic conicB(smajor*1.5, sminor, xcen+215, ycen, angle);
  // Eigen::Vector3d g; g.fill(0);
  // Eigen::Vector3d h; h.fill(0);
  std::tuple<Eigen::Vector3d, Eigen::Vector3d> gh;
  Eigen::Vector3d g_line, h_line;
  bool success = conicA.ConicIntersectionLines(conicB, gh);

  std::tie(g_line, h_line) = gh;

  std::vector<Conic> conics = {conicA, conicB};
  cv::Mat image(320, 480, CV_8UC3, 
                cv::Scalar(25, 25, 25));
  viz::drawEllipses(image, conics, viz::CV_colors);
  cv::Scalar red(0,0,255);
  viz::drawLine(image, h_line, cv::Scalar(0,0,255));
  // Showing image inside a window 
  cv::imshow("Conic Intersection", image); 
  // cv::waitKey(0); 
  ASSERT_TRUE(success);
}

TEST(ConicTest, CheckMatlab) {

  Conic conicA(245.848, 245.874, 1283.4, 1037.6, rad2deg(2.356));
  Conic conicB( 94.435, 261.000, 1808.5, 2081.3, rad2deg(1.120));
  // std::array<double,IMPLICIT_PARAM> implA = conicA.GetImplicit();
  // std::array<double,IMPLICIT_PARAM> implB = conicB.GetImplicit();
  // std::string iA = io::stringifyVector(implA, "Conic A: ");
  // std::string iB = io::stringifyVector(implB, "Conic B: ");
  // std::cout << iA << std::endl;
  // std::cout << iB << std::endl;


  std::tuple<Eigen::Vector3d, Eigen::Vector3d> gh;
  Eigen::Vector3d g_line, h_line;

  bool success = conicA.ConicIntersectionLines(conicB, gh);
  std::tie(g_line, h_line) = gh;
  // std::cout << "G line:\n" << g_line/g_line[1] << std::endl;
  // std::cout << "H line:\n" << h_line/h_line[1] << std::endl;
  ASSERT_TRUE(success);
}

TEST(ConicTest, ChooseCorrectIntersection) {
  Conic conicA(245.848, 245.874, 283.4, 0037.6, rad2deg(2.356));
  Conic conicB( 94.435, 261.000, 510.5, 0581.3, rad2deg(1.120));
  Eigen::Vector2d centerA, centerB;
  centerA = conicA.GetCenter();
  centerB = conicB.GetCenter();

  Eigen::Vector3d g, h;
  std::tuple<Eigen::Vector3d, Eigen::Vector3d> gh;
  bool success = conicA.ConicIntersectionLines(conicB, gh);
  EXPECT_TRUE(success);
  std::tie(g, h) = gh;
  // convert centers to homogeneous coordinates
  Eigen::Vector3d centerAHom = centerA.homogeneous();
  Eigen::Vector3d centerBHom = centerB.homogeneous();
  // get line connecting the two centers
  Eigen::Vector3d lineOfCenters = centerAHom.cross(centerBHom);
  // get point where lineOfCenters and g intersect
  Eigen::Vector2d gIntersect = g.cross(lineOfCenters).hnormalized();
  // get point where lineOfCenters and h intersect
  Eigen::Vector2d hIntersect = h.cross(lineOfCenters).hnormalized();

  double xmax, xmin, ymax, ymin;
  xmax = (centerA(0)>centerB(0)) ? centerA(0) : centerB(0);
  xmin = (centerA(0)<centerB(0)) ? centerA(0) : centerB(0);
  ymax = (centerA(1)>centerB(1)) ? centerA(1) : centerB(1);
  ymin = (centerA(1)<centerB(1)) ? centerA(1) : centerB(1);
  bool g_fits, h_fits;
  
  g_fits = gIntersect(0)>xmin && gIntersect(0)<xmax && gIntersect(1)>ymin && gIntersect(1)<ymax;
  h_fits = hIntersect(0)>xmin && hIntersect(0)<xmax && hIntersect(1)>ymin && hIntersect(1)<ymax;
  bool valid_intersection = false;
  if (g_fits ^ h_fits) {
    valid_intersection = true;
  }
  Eigen::Vector3d l_check;
  l_check = g_fits ? g : h;
  // std::cerr << "Result of intersection: " 
  //           << (g_fits ? "g":"h") << std::endl 
  //           << l_check << std::endl
  //           << "Wrong value:\n" << h << std::endl;
  EXPECT_TRUE(valid_intersection);

  Eigen::Vector3d l;
  bool match_success = ChooseIntersection(gh, centerA, centerB, l);
  EXPECT_TRUE(match_success);
  EXPECT_EQ(l_check, l);


  std::vector<Conic> conics = {conicA, conicB};
  cv::Mat image(720, 1280, CV_8UC3, 
                cv::Scalar(25, 25, 25));
  viz::drawEllipses(image, conics, viz::CV_colors);
  cv::Scalar   red(0,0,255);
  cv::Scalar green(0,255,0);
  viz::drawLine(image, l, green);
  viz::drawLine(image, h, red);
  // Showing image inside a window 
  cv::imshow("Choose Intersection", image); 
  // cv::waitKey(0); 
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

