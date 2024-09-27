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
    lat = 0, lon = 0, radius = 200;
    id = "defaultQuadric";
    quadric_default = new Quadric(lat, lon, radius, id);
  }

public:
  // ~QueueTest() override = default;
  double lat, lon, radius;
  std::string id;
  Quadric* quadric_default;
  Eigen::Vector3d origin = Eigen::Vector3d::Zero();
};

TEST_F(QuadricTest, MakingQuadric) {
  Eigen::Matrix3d t_e2m;
  // For clarity: we use a passive transformation, not an active rotation
  //   local x col vector is component in the Y-direction
  //   local y col vector is component in the Z-direction
  //   local z col vector is component in the X-direction
  t_e2m <<  0, 0, 1, 
            1, 0, 0, 
            0, 1, 0;
  ASSERT_TRUE(quadric_default->getQuadricTransformationMatrix().isApprox(t_e2m));
  ASSERT_EQ(quadric_default->getID(), id);
}

TEST_F(QuadricTest, InvalidCraterRadius) {
  double invalid_radius = 17000;
  double zero_radius = 0;
  std::string id_excess = "excessRadius";
  std::string id_zero = "zeroRadius";
  EXPECT_THROW(Quadric quad(lat, lon, invalid_radius, id_excess), std::runtime_error);
  EXPECT_THROW(Quadric quad(lat, lon, zero_radius, id_zero), std::runtime_error);
}

TEST_F(QuadricTest, InvalidPositionWithDependentNormal) {
  std::string id = "dependentNormal";
  EXPECT_THROW(Quadric quad(origin, radius, id), std::runtime_error);
}

TEST_F(QuadricTest, InvalidSurfaceNormal) {
  Eigen::Vector3d position = 1e5*latlon2bearing(lat, lon);
  Eigen::Vector3d surface_normal(0, 0, 0);
  std::string id = "zeroNormal";
  EXPECT_THROW(Quadric quad(position, radius, surface_normal, id), std::runtime_error);
}

TEST_F(QuadricTest, AngleBetweenQuadrics) {
  double lat1 =  0, lon1 = 0, radius1 = 50;
  double lat2 = 30, lon2 = 0, radius2 = 300;
  Quadric q1(lat1, lon1, radius1);
  Quadric q2(lat2, lon2, radius2);
  Quadric q1_copy(lat1, lon1, radius1);
  double angle12 = q1.getAngleBetweenQuadrics(q2);
  ASSERT_DOUBLE_EQ(rad2deg(angle12), std::abs(lat2-lat1));
  double angle_zero = q1.getAngleBetweenQuadrics(q1_copy);
  ASSERT_DOUBLE_EQ(angle_zero, 0.);
}

TEST_F(QuadricTest, AxisOfRotationQuadrics) {
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

TEST_F(QuadricTest, QuadricFrom3Points) {

  Eigen::Vector3d p1 = {10,  10, 0};
  Eigen::Vector3d p2 = {10,  00, 10};
  Eigen::Vector3d p3 = {10, -10, 0};

  Quadric quad(p1, p2, p3, "threePoints");
  Eigen::Vector3d center = quad.getLocation();
  double radius = quad.getRadius();
  ASSERT_DOUBLE_EQ(radius, 10);
  Eigen::Vector3d center_check = {10, 0, 0};
  ASSERT_TRUE(center.isApprox(center_check));
  Eigen::Vector3d normal_check = Eigen::Vector3d::UnitX();
  EXPECT_TRUE(quad.getPlane().normal().isApprox(normal_check) || 
              quad.getPlane().normal().isApprox(-normal_check));
}

TEST_F(QuadricTest, SamePlane) {
  Eigen::Vector3d p1 = {10,  10,   0};
  Eigen::Vector3d p2 = {10,   0,  10};
  Eigen::Vector3d p3 = {10, -10,   0};
  Eigen::Vector3d p4 = {10,   0, -10.00000001};

  Quadric quad1(p1, p2, p3, "threePoints1");
  Quadric quad2(p1, p2, p4, "threePoints2");
  ASSERT_TRUE(isSamePlane(quad1.getPlane(), quad2.getPlane()));
}

TEST_F(QuadricTest, RadiusCheckFromPoints) {
  double r = 100;
  double rho = calculateCraterRimFromRadius(r);
  double phi0 = deg2rad(10.), phi1 = deg2rad(115.), phi2 = deg2rad(280.);
  Eigen::Vector3d p1 = {rho, r*std::cos(phi0), r*std::sin(phi0)};
  Eigen::Vector3d p2 = {rho, r*std::cos(phi1), r*std::sin(phi1)};
  Eigen::Vector3d p3 = {rho, r*std::cos(phi2), r*std::sin(phi2)};

  Quadric quad(p1, p2, p3, "fromRadius");
  Quadric quad_check(0, 0, r, "pts_check");

  EXPECT_DOUBLE_EQ(r, quad.getRadius());

  EXPECT_TRUE(isSamePlane(quad, quad_check));
  EXPECT_TRUE(quad == quad_check);

}

TEST_F(QuadricTest, InPlaneXYQuadric) {
  double r = 1e4;
  double z = 0;
  double phi0 = deg2rad(10.), phi1 = deg2rad(115.), phi2 = deg2rad(280.);
  Eigen::Vector3d p1 = {r*std::sin(phi0), r*std::cos(phi0), z};
  Eigen::Vector3d p2 = {r*std::sin(phi1), r*std::cos(phi1), z};
  Eigen::Vector3d p3 = {r*std::sin(phi2), r*std::cos(phi2), z};
  Quadric quad_inplane(p1, p2, p3, "inPlane");
  Eigen::Hyperplane<double, 3> plane(Eigen::Vector3d(0,0,-1),0);
  EXPECT_TRUE(isSamePlane(plane, quad_inplane.getPlane()));
}

TEST_F(QuadricTest, InPlaneXZQuadric) {
  double r = 1e4;
  double y = 0;
  // TODO: there are some instances where points like (1,0,0),(0,1,0), and (-1,0,0)
  // fail to form a valid plane, which is strange
  Eigen::Vector3d p1 = {0, y,  r};
  Eigen::Vector3d p2 = {r, y,  0};
  Eigen::Vector3d p3 = {0, y, -r};
  Quadric quad_inplane(p1, p2, p3, "inStrangePlane");
  Eigen::Hyperplane<double, 3> plane(Eigen::Vector3d(0,-1,0),0);
  EXPECT_TRUE(isSamePlane(plane, quad_inplane.getPlane()));
}
