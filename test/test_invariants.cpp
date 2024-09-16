#include "gtest/gtest.h"
#include <vector>
#include <tuple>
#include <iostream>

#include "conics.h"
#include "io.h"
#include "opencv2/core/types.hpp"
#include "opencv2/imgproc.hpp"
#include "visuals.h"

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

TEST(InvariantTest, IntersectionLines) {
  Conic conicA(10, 7, 300, 50, 0);
  Conic conicB(15, 12, 100, 200, 0);
  Conic conicC(12, 8, 50, 200, 0);
  std::tuple<Eigen::Vector3d, Eigen::Vector3d> gh_ij, gh_jk, gh_ki;
  // Eigen::Vector3d lij, ljk, lki;
  // double invA, invB, invC;
  Eigen::Matrix3d locusA, locusB, locusC;
  // Eigen::Vector2d centerA, centerB, centerC;
  locusA = conicA.GetLocus();
  locusB = conicB.GetLocus();
  locusC = conicC.GetLocus();
  bool success_ab = IntersectionLines(locusA, locusB, gh_ij);
  EXPECT_TRUE(success_ab);
  bool success_bc = IntersectionLines(locusB, locusC, gh_jk);
  EXPECT_TRUE(success_bc);
  bool success_ca = IntersectionLines(locusC, locusA, gh_ki);
  EXPECT_TRUE(success_ca);
}

TEST(InvariantTest, InvariantTriad) {
  Conic conicA(10, 7, 300, 50, 0);
  Conic conicB(15, 12, 100, 200, 0);
  Conic conicC(12, 8, 50, 200, 0);
  Conic conicD(12, 8, 400, 20, 0);

  std::vector<Conic> conics = {conicA, conicB, conicC, conicD};

  cv::Mat image(500, 500, CV_8UC3, 
                cv::Scalar(25, 25, 25));
  viz::plotEllipses(image, conics, viz::CV_colors);
  Eigen::Vector3d my_line;
  my_line << 1.,-2.,3.;
  viz::plotline(image, my_line, cv::viz::Color::blue());
  // Showing image inside a window 
  cv::imshow("Output", image); 
  // cv::waitKey(0); 
  
  // Invariants
  std::array<double, NONCOPLANAR_INVARIANTS> invariantsABC, invariantsBCD, invariantsCDA, invariantsDAB;
  if(!computeCraterTriadInvariants(conicA, conicB, conicC, invariantsABC)) {
    std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
  }
  if(!computeCraterTriadInvariants(conicB, conicC, conicD, invariantsBCD)) {
    std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
  }
  if(!computeCraterTriadInvariants(conicC, conicD, conicA, invariantsCDA)) {
    std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
  }
  if(!computeCraterTriadInvariants(conicD, conicA, conicB, invariantsDAB)) {
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
}

