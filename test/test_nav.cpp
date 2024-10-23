#include "gtest/gtest.h"
#include <eigen3/Eigen/Geometry>
// #include <opencv2/imgproc.hpp> 
// #include <opencv2/highgui/highgui.hpp> 
// #include <random>

#include "camera.h"
#include "quadrics.h"
#include "conics.h"

class NavigationTest : public testing::Test {
protected:
  NavigationTest() {
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
    double dx, dy, skew, im_height, im_width, up, vp;
    cv::Size2i image_size;
    Eigen::Quaterniond quat = Eigen::Quaterniond::Identity();
    Eigen::Vector3d position{1,2,-3e4};

    Camera* cam;
    cv::viz::Camera* cv_cam;

};

Eigen::Matrix3d canonical(const Eigen::Matrix3d& image_conic) {
  double focal_length = 1;
  // Applying the concept of the canonical frame from Kanatani:
  // "3D interpretation of conics and orthogonality"
  double A = image_conic(0, 0);
  double B = image_conic(0, 1);
  double C = image_conic(1, 1);
  // D = image_conic[0, 2]/focal_length
  // E = image_conic[1, 2]/focal_length
  // F = image_conic[2, 2]/focal_length**2
  // if AC-B^2 > 0:
  //   - if A+C > 0, ellipse
  //   - if A+C < 0, imaginary conic
  // if AC-B^2  < 0, hyperbola
  // if AC-B^2 == 0, parabola
  double ACmB2 = A*C - std::pow(B,2);
  // this value must not be zero for an ellipse
  assert(std::abs(ACmB2)>1e-40);
  double ApC = A+C;
  assert( ApC > 0);
  double lambda1 = (ApC - std::sqrt(std::pow(ApC,2) - 4*ACmB2))/2;
  double lambda2 = (ApC + std::sqrt(std::pow(ApC,2) - 4*ACmB2))/2;
  double mu = std::pow(focal_length,2)/ACmB2;
  // to be an ellipse, mu*lambda_i must be > 0
  assert(mu*lambda1 > 0 and mu*lambda2 > 0);

  double a = std::sqrt(mu/lambda1);
  double b = std::sqrt(mu/lambda2);
  double a2 = std::pow(a,2);
  double b2 = std::pow(b,2);
  double A_prime = b2;
  double C_prime = a2;
  double F_prime = -a2*b2;
  Eigen::Matrix3d canonical = Eigen::Matrix3d::Identity();
  canonical.diagonal() = Eigen::Vector3d(A_prime, C_prime, F_prime);
  return canonical;
}

template<typename T>
std::vector<size_t> argsort(const std::vector<T> &array) {
    std::vector<size_t> indices(array.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&array](int left, int right) -> bool {
                  // sort indices according to corresponding array element
                  return array[left] < array[right];
              });

    return indices;
}

void conicBackprojection(const Eigen::Matrix3d& conic_locus, const double radius, Eigen::Vector3d& normal, double& dist) {
  // Applying the concept of the canonical frame from Kanatani:
  // "3D interpretation of conics and orthogonality"
  // Eigenvalue decomposition
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(conic_locus);
  if (eigensolver.info() != Eigen::Success) {
    std::cerr << "Eigenvalue decomposition failed!" << std::endl;
    return;
  }
  // std::cout << __func__ << "--> " << Conic(conic) << std::endl;
  
  // Eigenvalues and eigenvectors
  Eigen::Vector3d eigenvalues  = eigensolver.eigenvalues();
  Eigen::Matrix3d eigenvectors = eigensolver.eigenvectors();
  std::cout << "Eigenvalues: " << eigenvalues(0) << ", " << eigenvalues(1) << ", " << eigenvalues(2) << std::endl;
  std::cout << "Eigenvectors:\n" << eigenvectors << std::endl;
  
  // We (possibly) need to ensure that we have exactly two positive eigenvalues, else flip sign
  if((eigenvalues(0) * eigenvalues(1) * eigenvalues(2)) > 0) {
    eigenvalues *= -1.0;
    // std::swap(eigenvalues(0), eigenvalues(2));
    // eigenvectors.col(0).swap(eigenvectors.col(2));
  }
  std::vector<double> eigen_vec(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
  std::vector<size_t> indices = argsort(eigen_vec);
  std::cout << "******SORT INDICES: " << indices.at(0) << ", "<< indices.at(1) << ", "<< indices.at(2) << std::endl;
  // eigenvalues = eigenvalues(indices);
  Eigen::Vector3d sorted_eigenvalues = eigenvalues(indices);
  Eigen::Matrix3d sorted_eigenvectors;
  sorted_eigenvectors <<  eigenvectors.col(indices.at(0)), 
                          eigenvectors.col(indices.at(1)), 
                          eigenvectors.col(indices.at(2));
  std::cout << "-->Eigenvalues: " << sorted_eigenvalues.transpose() << std::endl;
  std::cout << "-->Eigenvectors:\n" << sorted_eigenvectors << std::endl;

  // lambda3 < 0 < lambda1 <= lambda2
  // u2 -> lambda2 | u3 -> lambda3
  Eigen::Vector3d::Index i1=1, i2=2, i3=0;
  double lambda1 = eigenvalues(i1);
  double lambda2 = eigenvalues(i2);
  double lambda3 = eigenvalues(i3);
  Eigen::Vector3d u1 = eigenvectors.col(i1);
  Eigen::Vector3d u2 = eigenvectors.col(i2);
  Eigen::Vector3d u3 = eigenvectors.col(i3);
  Eigen::Matrix3d swapped;
  swapped << u1, u2, u3;
  std::cout << "\nLambdas: " << lambda1 << ", " << lambda2 << ", " << lambda3 << std::endl;
  std::cout << "Associated u_i:\n" << swapped << std::endl;
  assert(lambda3 < 0);
  assert(lambda1 > 0);
  assert(lambda3 < lambda2);
  assert(lambda1 < lambda2);
  double scalar1 = std::sqrt((lambda2 - lambda1)/(lambda2 - lambda3));
  double scalar2 = std::sqrt((lambda1 - lambda3)/(lambda2 - lambda3));
  assert(!std::isnan(scalar1));
  assert(!std::isnan(scalar2));
  normal = scalar1*u2 + scalar2*u3;
  std::cout << "Quadric radius: " << radius << std::endl;
  dist = std::pow(std::abs(lambda1), 1.5) * radius;
}

int findOppositeSignedValue(const std::vector<double>& vec) {
  double sign_val = std::accumulate(
    begin(vec), end(vec), 1.0, std::multiplies<double>());
  auto it = std::find_if(vec.begin(), vec.end(), [sign_val](double num) {
    return (num * sign_val) > 0;
  });
  return it != vec.end() ? it - vec.begin() : -1;
}

void conicBackprojectionShiu(const Eigen::Matrix3d& conic, const double radius, Eigen::Vector3d& normal, double& dist) {
  // Applying the concept of the canonical frame from Shiu:
  // "3D loc. of circular and spherical features by monocular ... vision"
  // Eigenvalue decomposition
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(conic);
  if (eigensolver.info() != Eigen::Success) {
    std::cerr << "Eigenvalue decomposition failed!" << std::endl;
    return;
  }
  // std::cout << __func__ << "--> " << Conic(conic) << std::endl;
  
  // Eigenvalues and eigenvectors
  Eigen::Vector3d eigenvalues  = eigensolver.eigenvalues();
  Eigen::Matrix3d eigenvectors = eigensolver.eigenvectors();
  std::vector<int> indices {0, 1, 2};
  std::cout << "Eigenvalues: " << eigenvalues.transpose() << std::endl;

  // Find the eigenvalue with the different sign
  // Eqns 3.1.2.5 - 3.1.2.7
  Eigen::Vector3d f_d;
  double mu_d;
  int mu_d_idx;
  std::vector<double> eigen_vec(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
  std::vector<double> eigen_cpy = eigen_vec;
  mu_d_idx = findOppositeSignedValue(eigen_vec);
  // eigen_cpy.erase(mu_d_idx);
  if(mu_d_idx == -1) {
    std::cerr << "Eigenvalue of opposite sign not found: " << eigenvalues.transpose() << std::endl;
    return;
  }
  mu_d = eigenvalues(mu_d_idx);
  Eigen::MatrixXd eigenvectors_same = eigenvectors;
  eigen_cpy.erase(std::remove(eigen_cpy.begin(), eigen_cpy.end(), mu_d), eigen_cpy.end());
  indices.erase(std::remove(indices.begin(), indices.end(), mu_d_idx), indices.end());
  f_d = eigenvectors.col(mu_d_idx);

  // double lambda1 = eigenvalues(i1);
  // double lambda2 = eigenvalues(i2);
  double lambda3 = mu_d;
  double dotted = f_d.dot(Eigen::Vector3d::UnitZ());
  assert(dotted != 0);
  Eigen::Vector3d e3 = dotted > 0 ? f_d : -f_d;
  std::cout << "Unique value: " << lambda3 << " at index " << mu_d_idx << std::endl;
  std::cout << "Remaining eigenvalues: " << eigen_cpy.at(0) << ", " << eigen_cpy.at(1) << std::endl;

  // Get remaining lambdas
  // Eqns 3.1.2.8 - 3.1.2.11
  double omega1, omega2, lambda2, lambda1;
  Eigen::Vector3d g1, g2, e2, e1;
  int indexMax = indices.at(1);
  int indexMin = indices.at(0);
  omega1 = eigen_cpy.at(indexMax);
  omega2 = eigen_cpy.at(indexMin);
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
  double abs_l1 = std::abs(lambda1);
  double abs_l2 = std::abs(lambda2);
  double abs_l3 = std::abs(lambda3);
  double scalar = 1/std::abs(lambda2);
  double r = scalar * std::sqrt((abs_l1*abs_l3*(abs_l2 + abs_l3))/(abs_l1 + abs_l3));
  std::cout << "r: " << r << " | R/r: " << radius/r << std::endl;
}

TEST_F(NavigationTest, ConicBackprojection) {
  // Camera cam;
  Quadric quad(-15,0,100,"backprojection");
  Eigen::Vector3d location = Eigen::Vector3d::Zero();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  cam->resetCameraState();
  cam->moveX(1e0);
  cam->pointTo(location, up_vector);
  cam->moveRelative(Eigen::Vector3d(1e4, 0, 0));
  
  Eigen::Vector3d normal;
  double dist;
  std::cout << quad << std::endl;
  std::cout << cam << std::endl;
  std::cout << cam->getAttitude().inverse() << std::endl;
  std::cout << "Extrinsic matrix:\n" << cam->getExtrinsicMatrix() << std::endl;
  std::cout << "Projection matrix:\n" << cam->getProjectionMatrix() << std::endl;
  double dist1 = (quad.getLocation() - cam->getPosition()).norm();
  std::cout << "Distance should be roughly " << dist1 << "km\n";
  std::cout << "Value of lambda1 should be " << std::pow(dist1/quad.getRadius(), 2.0/3.0) << std::endl;
  // Eigen::Matrix3d conic_envelope = quad.projectToConicEnvelope(cam->getProjectionMatrix());
  Eigen::Matrix3d conic_locus = quad.projectToPlaneLocus(cam->getProjectionMatrix());
  // Eigen::Matrix3d plane_locus = cam->getImagePlaneLocus(conic_locus);
  conicBackprojectionShiu(conic_locus, quad.getRadius(), normal, dist);
  // conicBackprojection(conic_locus, quad.getRadius(), normal, dist);
  std::cout << "Conic normal: " << normal.transpose() << " | distance: " << dist << std::endl;
}

TEST_F(NavigationTest, KanataniEllipseCheck) {
  // double semimajor = 100, semiminor = 40, xc = 40, yc = 120, phi = deg2rad(10);
  // Conic conic(semimajor, semiminor, xc, yc, phi);
}

TEST_F(NavigationTest, ImageConic2PlaneConic) {
  Quadric quad0( 0,  0, 100, "Crater 1");
  // Quadric quad1( 10,  15, 50, "Crater 1");
  // Quadric quad2(-10,  15, 100, "Crater 2");
  // Quadric quad3(-30, -10, 75, "Crater 3");
  // Quadric quad4(  0, -20, 150, "Crater 4");

  Eigen::Vector3d up_axis(0,0,1);
  cam->move(Eigen::Vector3d(1e4, 0, 0));
  cam->pointTo(quad0.getLocation(), up_axis);
  

  Eigen::Matrix3d c_envelope = cam->projectQuadric(quad0.getEnvelope());
  Eigen::Matrix3d c_locus = adjugate(c_envelope);
  Conic image_conic(c_locus);
  Conic plane_conic(cam->getImagePlaneLocus(c_locus));
  Eigen::Vector2d ellipse_center = image_conic.getCenter();
  Eigen::Vector2d plane_center = plane_conic.getCenter();
  Eigen::Vector2d image_center = cam->getImageMidpoint();
  Eigen::Vector2d zeros = Eigen::Vector2d::Zero();
  ASSERT_TRUE(image_center.isApprox(ellipse_center, 1e-8));
  // TODO: For some reason the next line fails, but the subsequent line passes
  // ASSERT_TRUE(zeros.isApprox(plane_center, 1e-8));
  ASSERT_LE((plane_center-zeros).norm(), 1e-9);
}
