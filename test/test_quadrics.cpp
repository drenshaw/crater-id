#include "gtest/gtest.h"
#include <eigen3/Eigen/Geometry>
// #include <random>

#include "quadrics.h"
#include "io.h"
#include "math_utils.h"

// #include <iomanip>

/* // CATCH2 Testing Framework
#include <catch2/catch_test_macros.hpp>
#include <stdexcept>

TEST_CASE( "Factorial of 0 is 1 (fail)", "[single-file]" ) {
    REQUIRE( 1 == 1 );
    std::cout << "Test1\n";
}
*/

class QuadricTest : public testing::Test {
protected:
  QuadricTest() {
  }

public:
  // // ~QueueTest() override = default;
  // Queue<int> q0_;
  // Queue<int> q1_;
  // Queue<int> q2_;
};

TEST(QuadricTest, MakingQuadric) {
  double lat = 0, lon = 0, radius = 200;
  std::string id = "testLatLon";
  Quadric quad(lat, lon, radius, id);
  // std::cerr << quad << std::endl;
  // std::cerr << quad.GetLocus() << std::endl;
  Eigen::Matrix3d t_e2m;
  // For clarity: we use a passive transformation, not an active rotation
  //   local x col vector is component in the Y-direction
  //   local y col vector is component in the Z-direction
  //   local z col vector is component in the X-direction
  t_e2m <<  0, 0, 1, 
            1, 0, 0, 
            0, 1, 0;
  ASSERT_TRUE(quad.getQuadricTransformationMatrix()==t_e2m);
  // Eigen::Vector3d X = Eigen::Vector3d::UnitX();
  // Eigen::Vector3d Y = Eigen::Vector3d::UnitY();
  // Eigen::Vector3d Z = Eigen::Vector3d::UnitZ();
  // std::cout << "X vector pointing here : " << (t_e2m*X).transpose() << std::endl;
  // std::cout << "Y vector pointing here : " << (t_e2m*Y).transpose() << std::endl;
  // std::cout << "Z vector pointing here : " << (t_e2m*Z).transpose() << std::endl;
  
  Quadric north_pole(89, 0, radius, "northPole");
  // std::cerr << north_pole << std::endl;
  // std::cerr << north_pole.GetLocus() << std::endl;
  // std::cerr << north_pole.getQuadricTransformationMatrix() << std::endl;
  // std::cerr << "Center to crater radius: " << calculateCraterRimFromRadius(500) << std::endl;
}

TEST(QuadricTest, InvalidCraterRadius) {
  double lat = 0, lon = 0, invalid_radius = 17000;
  std::string id_excess = "excessRadius";
  EXPECT_THROW(Quadric quad(lat, lon, invalid_radius, id_excess), std::runtime_error);
  EXPECT_THROW(Quadric quad(lat, lon, 0, "zeroRadius"), std::runtime_error);
}

TEST(QuadricTest, InvalidPositionWithDependentNormal) {
  double radius = 100;
  Eigen::Vector3d position;
  position << 0, 0, 0;
  std::string id = "dependentNormal";
  EXPECT_THROW(Quadric quad(position, radius, id), std::runtime_error);
}

TEST(QuadricTest, InvalidSurfaceNormal) {
  double lat = 0, lon = 0, radius = 17000;
  Eigen::Vector3d position = 1e5*latlon2bearing(lat, lon);
  Eigen::Vector3d surface_normal;
  surface_normal << 0, 0, 0;
  std::string id = "zeroNormal";
  EXPECT_THROW(Quadric quad(position, radius, surface_normal, id), std::runtime_error);
}

TEST(QuadricTest, AngleBetweenQuadrics) {
  double lat1 =  0, lon1 = 0, radius1 = 50;
  double lat2 = 30, lon2 = 0, radius2 = 300;
  Quadric q1(lat1, lon1, radius1);
  Quadric q2(lat2, lon2, radius2);
  Quadric q1_copy(lat1, lon1, radius1);
  double angle12 = q1.getAngleBetweenQuadrics(q2);
  ASSERT_DOUBLE_EQ(rad2deg(angle12), std::abs(lat2-lat1));
  double angle_zero = q1.getAngleBetweenQuadrics(q1_copy);
  ASSERT_DOUBLE_EQ(angle_zero, 0.);
  Eigen::Hyperplane<double, 3> hplane1 = q1.getPlane();
  Eigen::Hyperplane<double, 3> hplane2 = q2.getPlane();
  

    using Line2 = Eigen::Hyperplane<float,2>;
    using Vec2  = Eigen::Vector2f;
    Vec2 a(8,2);
    Vec2 b(9,5);
    Vec2 c(6,6);
    Vec2 d(5,9);

    Line2 ac = Line2::Through(a,c);
    Line2 bd = Line2::Through(b,d);

    // std::cout << "Intersection:\n" << ac.intersection(bd) << '\n';
  // Eigen::Vector3d x = A.jacobiSvd(ComputeFullU|ComputeFullV).solve(b);
}

TEST(QuadricTest, AxisOfRotationQuadrics) {
  double lat1 =  0, lon1 = 0, radius1 = 50;
  double lat2 = 30, lon2 = 0, radius2 = 50;
  Quadric q1(lat1, lon1, radius1);
  Quadric q2(lat2, lon2, radius2);
  Quadric q1_copy(lat1, lon1, radius1);
  Eigen::Vector3d axis_normal = q1.getAxisNormalToQuadrics(q2);
  Eigen::Vector3d axis_check(3);
  // Moving from 0 degrees latitude up to 30 degrees latitude makes axis in -y direction
  axis_check << 0, -1, 0;
  ASSERT_TRUE(axis_normal.isApprox(axis_check));
  ASSERT_THROW(q1.getAxisNormalToQuadrics(q1_copy), std::runtime_error);
}

TEST(QuadricTest, QuadricFrom3Points) {

  Eigen::Vector3d p1 = {10,  10, 0};
  Eigen::Vector3d p2 = {10,  00, 10};
  Eigen::Vector3d p3 = {10, -10, 0};

  Quadric quad(p1, p2, p3, "threePoints");
  std::cout << quad << std::endl;
  Eigen::Vector3d center = quad.getLocation();
  double radius = quad.getRadius();
  ASSERT_DOUBLE_EQ(radius, 10);
  Eigen::Vector3d center_check = {10, 0, 0};
  ASSERT_TRUE(center.isApprox(center_check));
  Eigen::Vector3d normal_check = Eigen::Vector3d::UnitX();
  EXPECT_TRUE(quad.getPlane().normal().isApprox(normal_check) || 
              quad.getPlane().normal().isApprox(-normal_check));
  // std::cout << "Locus:\n" << quad.getLocus() << std::endl;
}

TEST(QuadricTest, RadiusCheckFromPoints) {
  double radius = 100;
  // double rho = std::sqrt(std::pow(R_MOON, 2) - std::pow(radius, 2));
  double z = R_MOON * (M_PI_2 - std::asin(radius/R_MOON));

  Eigen::Vector3d p1 = { radius,      0, z};
  Eigen::Vector3d p2 = {     0,  radius, z};
  Eigen::Vector3d p3 = {-radius,      0, z};

  Quadric quad(p1, p2, p3, "fromRadius");
  EXPECT_EQ(radius, quad.getRadius());
  Eigen::Vector3d normal_check = Eigen::Vector3d::UnitZ();
  EXPECT_TRUE(quad.getPlane().normal().isApprox(normal_check) || 
              quad.getPlane().normal().isApprox(-normal_check));

}
