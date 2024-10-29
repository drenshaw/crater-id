#include "gtest/gtest.h"
#include <eigen3/Eigen/Geometry>
// #include <opencv2/imgproc.hpp> 
// #include <opencv2/highgui/highgui.hpp> 
// #include <random>

#include "navigation.h"
#include "camera.h"
#include "quadrics.h"
#include "conics.h"

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
  Eigen::Matrix3d Kinv = cam->getInverseIntrinsicMatrix();
  
  std::array<Eigen::Vector3d, 2> centers, normals;
  Eigen::Vector3d center1, center2, normal1, normal2;
  std::cout << quadric_default << std::endl;
  std::cout << cam << std::endl;

  Eigen::Matrix3d c_locus = cam->projectQuadricToLocus(quadric_default->getLocus());
  Eigen::Matrix3d plane_locus = Kinv * c_locus * Kinv.transpose();
  // double dist;
  // double dist2center = (quad.getLocation() - cam->getPosition()).norm();
  // Eigen::Matrix3d conic_locus = quad.projectToPlaneLocus(extrinsic);
  Eigen::Matrix3d conic_locus = quadric_default->projectToPlaneLocus(proj);
  std::cout << "Plane:\n" << quadric_default->getLocus() << std::endl;
  std::cout << "Plane1:\n" << conic_locus << std::endl;
  Conic cc(conic_locus);
  std::cout << cc << std::endl;
  Shiu::conicBackprojection(conic_locus, quadric_default->getRadius(), centers, normals);
  // conicBackprojection(conic_locus, quad.getRadius(), normal, dist);
  std::cout << "Conic center1: " << centers.at(0).transpose() << " | normal1: " << normals.at(0).transpose() << std::endl;
  std::cout << "Conic center2: " << centers.at(1).transpose() << " | normal2: " << normals.at(1).transpose() << std::endl;
}

void getBackprojectionRemainingLambdas( const std::vector<double>& eigenvalues,
                                            const Eigen::Matrix3d& eigenvectors,
                                            const Eigen::Vector3d e3,
                                            double& lambda1, double& lambda2,
                                            Eigen::Vector3d& e1, Eigen::Vector3d& e2) {
  double omega1, omega2;
  Eigen::Vector3d g1, g2;
  // Get remaining lambdas
  // Eqns 3.1.2.8 - 3.1.2.11
  
  Eigen::Vector2d eig;
  eig << eigenvalues.at(0), eigenvalues.at(1);
  Eigen::Vector2d::Index maxRow;
  eig.maxCoeff(&maxRow);
  std::cout << "Max row: " << maxRow << std::endl;
  int indexMin = 0;
  int indexMax = 1;
  omega1 = eigenvalues.at(indexMax);
  omega2 = eigenvalues.at(indexMin);
  g1 = eigenvectors.col(indexMax);
  g2 = eigenvectors.col(indexMin);
  if(omega1 == omega2) {
    std::cerr << "Not confident in handling right cones yet!\n";
  }
  if(std::abs(omega1) < std::abs(omega2)) {
    lambda2 = omega1;
    lambda1 = omega2;
    e2 = g1;
  }
  else {
    lambda2 = omega2;
    lambda1 = omega1;
    e2 = g2;
  }
  e1 = e2.cross(e3);
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
  
  // Store eigenvalues in a vector for easy manipulation
  std::vector<double> eigenvalues(eigenval.data(), eigenval.data() + eigenval.size());

  double lambda3, lambda2, lambda1;
  Eigen::Vector3d g1, g2, u3, u2, u1;
  int mu_d_idx;  

  // Find the eigenvalue with the different sign
  // Eqns 3.1.2.5 - 3.1.2.7
  Shiu::getBackprojectionLambda1(eigenvalues, eigenvec, mu_d_idx, u3);
  lambda3 = eigenval(mu_d_idx);
  eigenvalues.erase(eigenvalues.begin() + mu_d_idx);

  // std::cout << "Unique value: " << lambda3 << " at index " << mu_d_idx << std::endl;
  assert(eigenvalues.size() == 2);
  // std::cout << "Remaining eigenvalues: " << eigenvalues.at(0) << ", " << eigenvalues.at(1) << std::endl;
  assert(eigenvalues.at(0) * eigenvalues.at(1) > 0);

  // Get remaining lambdas
  // Eqns 3.1.2.8 - 3.1.2.11
  getBackprojectionRemainingLambdas(eigenvalues, eigenvec, u3, lambda1, lambda2, u1, u2);
}

TEST_F(NavigationTest, ChristianBackprojection) {
  // Camera cam;
  // Quadric quad(-30,0,100,"backprojection");
  // Eigen::Vector3d location = Eigen::Vector3d::Zero();
  Eigen::Vector3d location = R_MOON*Eigen::Vector3d::UnitZ();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  // cam->resetCameraState();
  cam->moveX(1e4);
  cam->pointTo(location, up_vector);
  
  // Eigen::MatrixXd proj = cam->getProjectionMatrix();
  Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();
  
  std::array<Eigen::Vector3d, 2> centers, normals;
  Eigen::Vector3d center1, center2, normal1, normal2;
  // double dist;
  // double dist2center = (quad.getLocation() - cam->getPosition()).norm();
  // Eigen::Matrix3d conic_locus = quad.projectToPlaneLocus(extrinsic);
  Eigen::Matrix3d conic_locus = quadric_default->projectToPlaneLocus(extrinsic);
  conicBackprojection(conic_locus, quadric_default->getRadius(), centers, normals);
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
