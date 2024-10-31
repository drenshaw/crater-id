#include "gtest/gtest.h"
#include <eigen3/Eigen/Geometry>
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
// #include <random>

#include "navigation.h"
#include "camera.h"
#include "quadrics.h"
#include "conics.h"
#include "visuals.h"

class NavigationTest : public testing::Test {
protected:
  NavigationTest() {
    latitude = -30, longitude = 0, radius = 200;
    id = "defaultQuadric";
    quadric_default = new Quadric(latitude, longitude, radius, id);
    dx = 1000, dy = 1000, skew = 0, im_height = 2048, im_width = 2592;
    up = (im_width+1)/2;
    vp = (im_height+1)/2;
    image_size = cv::Size(im_width, im_height);
    cam = new Camera(dx, dy, up, vp, skew, image_size, quat, position);
    cv_cam = new cv::viz::Camera(dx, dy, up, vp, image_size);
  }
  ~NavigationTest() override {
    // delete quadric_default;
      delete cam;
      delete cv_cam;
  }
  public:
  double latitude, longitude, radius;
  std::string id;
  Quadric* quadric_default;
    double dx, dy, skew, im_height, im_width, up, vp;
    cv::Size2i image_size;
    Eigen::Quaterniond quat = Eigen::Quaterniond::Identity();
    Eigen::Vector3d position{1e4, 0, 0};

    Camera* cam;
    cv::viz::Camera* cv_cam;

};

TEST_F(NavigationTest, ShiuConicBackprojection) {
  // Camera cam;
  // Quadric quad(-30,0,100,"backprojection");
  // Eigen::Vector3d location = Eigen::Vector3d::Zero();
  Eigen::Vector3d location = R_MOON*Eigen::Vector3d::UnitZ();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  // cam->resetCameraState();
  cam->moveX(1e4);
  cam->pointTo(location, up_vector);
  
  Eigen::MatrixXd proj = cam->getProjectionMatrix();
  Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();
  // Eigen::Matrix3d Kinv = cam->getInverseIntrinsicMatrix();
  
  std::array<Eigen::Vector3d, 2> centers, normals;
  Eigen::Vector3d center1, center2, normal1, normal2;
  // std::cout << *quadric_default << std::endl;
  // std::cout << cam << std::endl;

  // Eigen::Matrix3d c_locus = cam->projectQuadricToLocus(quadric_default->getLocus());
  // Eigen::Matrix3d plane_locus = Kinv * c_locus * Kinv.transpose();
  // double dist;
  // double dist2center = (quad.getLocation() - cam->getPosition()).norm();
  // Eigen::Matrix3d conic_locus = quad.projectToPlaneLocus(extrinsic);
  Eigen::Matrix3d conic_locus = quadric_default->projectToPlaneLocus(extrinsic);
  // std::cout << "Plane:\n" << quadric_default->getLocus() << std::endl;
  // std::cout << "Plane1:\n" << conic_locus << std::endl;
  Conic cc(conic_locus);
  std::cout << cc << std::endl;
  Shiu::conicBackprojection(conic_locus, quadric_default->getRadius(), centers, normals);
}

TEST_F(NavigationTest, ChristianBackprojection) {
  // Camera cam;
  // Quadric quad(-30,0,100,"backprojection");
  // Eigen::Vector3d location = Eigen::Vector3d::Zero();
  Eigen::Vector3d look_here = -1e2*Eigen::Vector3d::UnitZ();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  // cam->resetCameraState();
  cam->moveX(1e4);
  cam->moveY(-1e3);
  cam->moveZ(-2e3);
  cam->pointTo(look_here, up_vector);
  
  // Eigen::MatrixXd proj = cam->getProjectionMatrix();
  Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();

  Eigen::Matrix3d att = cam->getAttitudeMatrix();
  Eigen::Vector3d centerWrtCamWorld = quadric_default->getLocation() - cam->getPosition();
  Eigen::Vector3d centerWrtCam = att * centerWrtCamWorld;
  double dist2center = (centerWrtCamWorld).norm();
  Eigen::Vector3d normal_test = att * quadric_default->getNormal();
  std::cout << "---\n---Truth Data---\n";
  std::cout << "---DISTANCE TO QUADRIC CENTER: " << dist2center << std::endl;
  std::cout << "---CAMERA TO QUADRIC CENTER: " << centerWrtCam.transpose() << std::endl;
  std::cout << "---QUADRIC NORMAL WRT CAMERA: " << normal_test.transpose() << std::endl;
  
  std::array<Eigen::Vector3d, 2> centers, normals;
  Eigen::Vector3d center1, center2, normal1, normal2;
  // double dist;
  // Eigen::Matrix3d conic_locus = quad.projectToPlaneLocus(extrinsic);
  Eigen::Matrix3d conic_locus = quadric_default->projectToPlaneLocus(extrinsic);
  conicBackprojection(conic_locus, quadric_default->getRadius(), centers, normals);
  Eigen::Vector3d q_normal = Eigen::Vector3d::UnitZ();
  center1 = centers.at(0);
  center2 = centers.at(1);
  normal1 = normals.at(0);
  normal2 = normals.at(1);
  std::cout << "+++\n+++Calculations+++\n";
  std::cout << "Center (cam) 1: " << center1.transpose() << std::endl;
  std::cout << "Center (cam) 2: " << center2.transpose() << std::endl;
  std::cout << "Normal (cam) 1: " << normal1.transpose() << std::endl;
  std::cout << "Normal (cam) 2: " << normal2.transpose() << std::endl;
  std::cout << "Angle 1: " << rad2deg(std::acos(q_normal.dot(-normal1))) << std::endl;
  std::cout << "Angle 2: " << rad2deg(std::acos(q_normal.dot(-normal2))) << std::endl;
  ASSERT_TRUE(centerWrtCam.isApprox(center1, 1e-3) || centerWrtCam.isApprox(center2, 1e-3));
  ASSERT_TRUE( normal_test.isApprox(normal1, 1e-3) ||  normal_test.isApprox(normal2, 1e-3));
  // conicBackprojection(conic_locus, quad.getRadius(), normal, dist);
}

TEST_F(NavigationTest, KanataniEllipseCheck) {
  // double semimajor = 100, semiminor = 40, xc = 40, yc = 120, phi = deg2rad(10);
  // Conic conic(semimajor, semiminor, xc, yc, phi);
}

TEST_F(NavigationTest, ImageConic2PlaneConic) {
  // Quadric quad0( 0,  0, 100, "Crater 1");
  // // Quadric quad1( 10,  15, 50, "Crater 1");
  // // Quadric quad2(-10,  15, 100, "Crater 2");
  // // Quadric quad3(-30, -10, 75, "Crater 3");
  // // Quadric quad4(  0, -20, 150, "Crater 4");

  // Eigen::Vector3d up_axis(0,0,1);
  // cam->move(Eigen::Vector3d(1e4, 0, 0));
  // cam->pointTo(quad0.getLocation(), up_axis);
  

  // Eigen::Matrix3d c_envelope = cam->projectQuadric(quad0.getEnvelope());
  // Eigen::Matrix3d c_locus = adjugate(c_envelope);
  // Conic image_conic(c_locus);
  // Conic plane_conic(cam->getImagePlaneLocus(c_locus));
  // Eigen::Vector2d ellipse_center = image_conic.getCenter();
  // Eigen::Vector2d plane_center = plane_conic.getCenter();
  // Eigen::Vector2d image_center = cam->getImageMidpoint();
  // Eigen::Vector2d zeros = Eigen::Vector2d::Zero();
  // ASSERT_TRUE(image_center.isApprox(ellipse_center, 1e-8));
  // // TODO: For some reason the next line fails, but the subsequent line passes
  // // ASSERT_TRUE(zeros.isApprox(plane_center, 1e-8));
  // ASSERT_LE((plane_center-zeros).norm(), 1e-9);
}
