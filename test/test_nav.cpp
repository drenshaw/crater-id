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

TEST_F(NavigationTest, ConicBackprojection) {
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
  // conicBackprojection(conic_locus, quad.getRadius(), normal, dist);
  std::cout << "Conic center1: " << centers.at(0).transpose() << " | normal1: " << normals.at(0).transpose() << std::endl;
  std::cout << "Conic center2: " << centers.at(1).transpose() << " | normal2: " << normals.at(1).transpose() << std::endl;
}

void conicBackprojection( const Eigen::Matrix3d& conic, const double radius, 
                          std::array<Eigen::Vector3d, 2>& centers, 
                          std::array<Eigen::Vector3d, 2>& normals) {
  // Applying the transforms found in Christian white paper
  // "Perspective projection of ellipses and ellipsoids with applications to spacecraft navigation"
  
  // Eigenvalues and eigenvectors
  Eigen::Vector3d eigenval;
  Eigen::Matrix3d eigenvec;
  if(!getEigenstuffConic(conic, eigenval, eigenvec)) {
    std::cerr << __func__ << "-> backprojection eigenstuff failed.\n";
    return;
  }

  double lambda3, lambda2, lambda1;
  Eigen::Vector3d g1, g2, u3, u2, u1;
  int mu_d_idx;  

  // Find the eigenvalue with the different sign
  // Eqns 3.1.2.5 - 3.1.2.7
  std::vector<int> indices = {0, 1, 2};
  Shiu::getBackprojectionLambda1(eigenval, eigenvec, mu_d_idx, u3);
  lambda3 = eigenval(mu_d_idx);
  // std::cout << "Eigenvalues: " << eigenval.transpose() << std::endl;
  // std::cout << "Eigenvectors:\n" << eigenvec << std::endl;
  indices.erase(indices.begin() + mu_d_idx);
  Eigen::Vector2d rem_eigenval = eigenval(indices);
  Eigen::MatrixXd rem_eigenvec = eigenvec(Eigen::all, indices);

  // Get remaining lambdas
  // Eqns 3.1.2.8 - 3.1.2.11
  Shiu::getBackprojectionLambda2(rem_eigenval, rem_eigenvec, lambda1, lambda2, u2);
  u1 = u2.cross(u3);
  Eigen::Matrix3d diagg;
  diagg << u1, u2, u3;
  // std::cout << "Diagonalized:\n" << diagg.transpose() * conic * diagg << std::endl;
  // X-axis corresponds to semiminor (smaller) axis
  double kx2 = 1/std::abs(lambda1);
  double ky2 = 1/std::abs(lambda2);
  double kz2 = 1/std::abs(lambda3);
  double alpha = (ky2 - kx2) / (ky2 + kz2);
  // double z_prime = 1/(1+std::sqrt(alpha));
  double kxz2 = kx2/kz2;
  double kxz = std::sqrt(kxz2);
  // double Xi = kxz/(1-alpha);
  // double Zi = 1/(1-alpha);
  double sqrtAlphaPlusKxz = std::sqrt(alpha + kxz2);
  // double rho = sqrtAlphaPlusKxz/(1 - alpha);
  Eigen::Vector3d center_posX, center_negX, normal_posX, normal_negX;
  Eigen::Vector3d center_pos, center_neg, normal_pos, normal_neg;
  double scalarR = radius/sqrtAlphaPlusKxz;
  double scalarN = 1/sqrtAlphaPlusKxz;
  center_posX <<  std::sqrt(alpha*kxz2), 0, 1;
  center_negX << -std::sqrt(alpha*kxz2), 0, 1;
  normal_posX <<  std::sqrt(alpha), 0, -kxz;
  normal_negX << -std::sqrt(alpha), 0, -kxz;
  center_posX *= scalarR;
  center_negX *= scalarR;
  normal_posX *= scalarN;
  normal_negX *= scalarN;
  center_pos  = diagg * center_posX;
  center_neg  = diagg * center_negX;
  normal_pos  = diagg * normal_posX;
  normal_neg  = diagg * normal_negX;
  centers = {center_pos, center_neg};
  normals = {normal_pos, normal_neg};
  // double dist = radius / rho;
  // std::cout << "XYZ: " << center_posX.transpose() << " | Normal: " << normal_posX.transpose() << " | dist: " << dist << std::endl;

  // // Eqn 3.1.3.6
  // // double dist = getBackprojectionDistance(lambda1, lambda2, lambda3, radius);
  // Shiu::getBackprojectionCenter(radius, lambda1, lambda2, lambda3, u1, u3, centers);
  // Shiu::getBackprojectionNormal(lambda1, lambda2, lambda3, u1, u3, normals);
  // Eigen::Vector3d center1, center2, normal1, normal2;
  // center1 = centers.at(0);
  // center2 = centers.at(1);
  // normal1 = normals.at(0);
  // normal2 = normals.at(1);
  // std::cout << "Center1: " << center1.transpose() << " | Distance: " << center1.norm() << std::endl;
  // std::cout << "Normal1: " << normal1.transpose() << std::endl;
  // std::cout << "Center2: " << center2.transpose() << " | Distance: " << center2.norm() << std::endl;
  // std::cout << "Normal2: " << normal2.transpose() << std::endl;
}

TEST_F(NavigationTest, ChristianBackprojection) {
  // Camera cam;
  // Quadric quad(-30,0,100,"backprojection");
  // Eigen::Vector3d location = Eigen::Vector3d::Zero();
  Eigen::Vector3d look_here = -1e3*Eigen::Vector3d::UnitZ();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  // cam->resetCameraState();
  cam->moveX(1e4);
  cam->pointTo(look_here, up_vector);
  
  // Eigen::MatrixXd proj = cam->getProjectionMatrix();
  Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();

  Eigen::Matrix3d att = cam->getAttitudeMatrix();
  Eigen::Vector3d centerWrtCamWorld = quadric_default->getLocation() - cam->getPosition();
  Eigen::Vector3d centerWrtCam = att * centerWrtCamWorld;
  double dist2center = (centerWrtCamWorld).norm();
  std::cout << "---\n---Truth Data---\n";
  std::cout << "---DISTANCE TO QUADRIC CENTER: " << dist2center << std::endl;
  std::cout << "---CAMERA TO QUADRIC CENTER: " << centerWrtCam.transpose() << std::endl;
  std::cout << "---QUADRIC NORMAL WRT CAMERA: " << (att * quadric_default->getNormal()).transpose() << std::endl;
  
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
