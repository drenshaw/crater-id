#include "io.h"
#include "conics.h"
#include "visuals.h"

#include "gtest/gtest.h"
#include <eigen3/Eigen/Dense>
#include <vector>
#include <tuple>
#include <iostream>

#include <opencv2/core/core.hpp> 
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include "opencv2/viz/types.hpp"

class InvariantTest : public testing::Test {
  protected:
    InvariantTest() {
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

// TEST(InvariantTest, IntersectionLines) {
//   Conic conicA(10, 7, 300, 50, 0);
//   Conic conicB(15, 12, 100, 200, 0);
//   Conic conicC(12, 8, 50, 200, 0);
//   std::tuple<Eigen::Vector3d, Eigen::Vector3d> gh_ij, gh_jk, gh_ki;
//   // Eigen::Vector3d lij, ljk, lki;
//   // double invA, invB, invC;
//   Eigen::Matrix3d locusA, locusB, locusC;
//   // Eigen::Vector2d centerA, centerB, centerC;
//   locusA = conicA.GetLocus();
//   locusB = conicB.GetLocus();
//   locusC = conicC.GetLocus();
//   bool success_ab = IntersectionLines(locusA, locusB, gh_ij);
//   EXPECT_TRUE(success_ab);
//   bool success_bc = IntersectionLines(locusB, locusC, gh_jk);
//   EXPECT_TRUE(success_bc);
//   bool success_ca = IntersectionLines(locusC, locusA, gh_ki);
//   EXPECT_TRUE(success_ca);

//   std::vector<Conic> conics = {conicA, conicB, conicC};
//   cv::Mat image(240, 320, CV_8UC3, 
//                 cv::Scalar(150, 150, 150));
//   std::vector<Eigen::Vector3d> lines;
//   lines = {gh_ij, gh_jk, gh_ki};
//   cv::Mat outImg;
//   viz::drawEllipses(image, conics, viz::CV_colors);
//   viz::drawLines(image, lines, viz::CV_colors);
//   cv::resize(image, outImg, cv::Size(), 0.25, 0.25);
//   cv::imshow("Invariants", outImg); 
//   cv::waitKey(0); 
// }

TEST(InvariantTest, InvariantTriad) {
  Conic conicA(245.848, 245.874, 1283.4, 1037.6, rad2deg(2.356));
  Conic conicB( 94.435, 261.000, 1808.5, 2081.3, rad2deg(1.120));
  Conic conicC(132.956, 210.487,  849.7, 1901.0, rad2deg(2.042));
  Conic conicD(132.956, 210.487, 1849.7, 0901.0, rad2deg(2.042));

  Eigen::Vector3d lij, lik, lil, ljk, ljl, lkl;
  if(!conicA.ChooseConicIntersection(conicB, lij)) {
    std::cerr<<"IntersectionLines error ij\n";
  }
  if(!conicA.ChooseConicIntersection(conicC, lik)) {
    std::cerr<<"IntersectionLines error jk\n";
  }
  if(!conicA.ChooseConicIntersection(conicD, lil)) {
    std::cerr<<"IntersectionLines error ki\n";
  }
  if(!conicB.ChooseConicIntersection(conicC, ljk)) {
    std::cerr<<"IntersectionLines error ki\n";
  }
  if(!conicB.ChooseConicIntersection(conicD, ljl)) {
    std::cerr<<"IntersectionLines error ki\n";
  }
  if(!conicC.ChooseConicIntersection(conicD, lkl)) {
    std::cerr<<"IntersectionLines error ki\n";
  }
  // Invariants
  double invA, invB, invC, invD;
  if(!invariants::computeInvariant(lij, lik, conicA.GetLocus(), invA)) {
    std::cerr<<"computeInvariant error A\n";
    ASSERT_TRUE(false);
  }
  if(!invariants::computeInvariant(ljk, lij, conicB.GetLocus(), invB)) {
    std::cerr<<"computeInvariant error B\n";
    ASSERT_TRUE(false);
  }
  if(!invariants::computeInvariant(lik, ljk, conicC.GetLocus(), invC)) {
    std::cerr<<"computeInvariant error C\n";
    ASSERT_TRUE(false);
  }
  if(!invariants::computeInvariant(lil, lkl, conicD.GetLocus(), invD)) {
    std::cerr<<"computeInvariant error C\n";
    ASSERT_TRUE(false);
  }
  double invA_check = 0.3556067233739398;
  double invB_check = 0.4065675587075122;
  double invC_check = 0.5830085148020281;


  EXPECT_NEAR(invA, invA_check, 1e-3);
  EXPECT_NEAR(invB, invB_check, 1e-3);
  EXPECT_NEAR(invC, invC_check, 1e-3);




  cv::Mat image(2400, 3200, CV_8UC3, 
                cv::Scalar(150, 150, 150));
  std::vector<Conic> conics = {
    conicA, 
    conicB, 
    conicC, 
    conicD
  };
  viz::drawEllipses(image, conics, viz::CV_colors);
  std::vector<Eigen::Vector3d> my_lines = {
    lij, 
    lik, 
    lil, 
    ljk, 
    ljl, 
    lkl
  };
  std::vector<std::string> my_text = {
    "lij", 
    "lik", 
    "lil", 
    "ljk", 
    "ljl", 
    "lkl"
  };
  viz::drawLines(image, my_lines, my_text, viz::CV_colors);
  // Showing image inside a window 
  cv::Mat outImg;
  cv::resize(image, outImg, cv::Size(), 0.4, 0.4);
  cv::imshow("Invariant Triad", outImg); 
  cv::waitKey(0); 
  
  // Invariants
  std::array<double, NONCOPLANAR_INVARIANTS> invariantsABC, invariantsBCD, invariantsCDA, invariantsDAB;
  if(!invariants::computeCraterTriadInvariants(conicA, conicB, conicC, invariantsABC)) {
    std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
  }
  if(!invariants::computeCraterTriadInvariants(conicB, conicC, conicD, invariantsBCD)) {
    std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
  }
  if(!invariants::computeCraterTriadInvariants(conicC, conicD, conicA, invariantsCDA)) {
    std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
  }
  if(!invariants::computeCraterTriadInvariants(conicD, conicA, conicB, invariantsDAB)) {
    std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
  }
  std::string invABC = io::stringifyVector(invariantsABC, "Invariants ABC: ");
  std::string invBCD = io::stringifyVector(invariantsBCD, "Invariants BCD: ");
  std::string invCDA = io::stringifyVector(invariantsCDA, "Invariants CDA: ");
  std::string invDAB = io::stringifyVector(invariantsDAB, "Invariants DAB: ");
  std::cout << invABC << std::endl;
  std::cout << invBCD << std::endl;
  std::cout << invCDA << std::endl;
  std::cout << invDAB << std::endl;
  if(!invariants::computeCraterTriadInvariants(conicA, conicB, conicC, invariantsABC)) {
    std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
  }
  EXPECT_DOUBLE_EQ(invA, invariantsABC.at(0));
  EXPECT_DOUBLE_EQ(invB, invariantsABC.at(1));
  EXPECT_DOUBLE_EQ(invC, invariantsABC.at(2));
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
  cv::Mat image(512, 640, CV_8UC3, 
                cv::Scalar(50, 50, 50));
  viz::drawEllipses(image, conics, viz::CV_colors);
  cv::Scalar red  = cv::viz::Color::red();
  cv::Scalar blue = cv::viz::Color::blue();
  viz::drawLine(image, g_line, "g_line", red);
  viz::drawLine(image, h_line, "h_line", blue);
  // Showing image inside a window 
  cv::imshow("Conic Intersection: `h` is the correct line", image); 
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
  Conic conicC(132.956, 210.487, 849.7, 1901.0, rad2deg(2.042));
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
  bool match_success = invariants::ChooseIntersection(gh, centerA, centerB, l);
  EXPECT_TRUE(match_success);
  EXPECT_EQ(l_check, l);


  std::vector<Conic> conics = {conicA, conicB};
  cv::Mat image(720, 1280, CV_8UC3, 
                cv::Scalar(25, 25, 25));
  viz::drawEllipses(image, conics, viz::CV_colors);
  cv::Scalar   red(0,0,255);
  cv::Scalar green(0,255,0);
  viz::drawLine(image, l, "true",  green);
  viz::drawLine(image, h, "false", red);
  // Showing image inside a window 
  cv::imshow("Choose Intersection", image); 
  // cv::waitKey(0); 
}

