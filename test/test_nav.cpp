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
  }
  ~NavigationTest() override {
    // delete quadric_default;
  }

public:
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

void conicBackprojection(const Eigen::Matrix3d& conic, const double radius, Eigen::Vector3d& normal, double& dist) {
  // Applying the concept of the canonical frame from Kanatani:
  // "3D interpretation of conics and orthogonality"
  // Eigenvalue decomposition
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(conic);
  if (eigensolver.info() != Eigen::Success) {
    std::cerr << "Eigenvalue decomposition failed!" << std::endl;
    return;
  }
  std::cout << __func__ << "--> " << Conic(conic) << std::endl;
  
  // Eigenvalues and eigenvectors
  Eigen::Vector3d eigenvalues = eigensolver.eigenvalues();
  Eigen::Matrix3d eigenvectors = eigensolver.eigenvectors();
  // lambda3 < 0 < lambda1 <= lambda2
  // u2 -> lambda2 | u3 -> lambda3
  Eigen::Vector3d::Index i1=1, i2=2, i3=0;
  double lambda1 = eigenvalues(i1);
  double lambda2 = eigenvalues(i2);
  double lambda3 = eigenvalues(i3);
  std::cout << "Eigenvalues: " << lambda1 << ", " << lambda2 << ", " << lambda3 << std::endl;
  std::cout << "Eigenvectors:\n" << eigenvectors << std::endl;
  assert(lambda3 < 0);
  assert(lambda1 > 0);
  assert(lambda3 < lambda2);
  assert(lambda1 < lambda2);
  Eigen::Vector3d u2 = eigenvectors.col(i2);
  Eigen::Vector3d u3 = eigenvectors.col(i3);
  double scalar1 = std::sqrt((lambda2 - lambda1)/(lambda2 - lambda3));
  double scalar2 = std::sqrt((lambda1 - lambda3)/(lambda2 - lambda3));
  assert(!std::isnan(scalar1));
  assert(!std::isnan(scalar2));
  normal = scalar1*u2 + scalar2*u3;
  std::cout << "Quadric radius: " << radius << std::endl;
  dist = std::pow(std::abs(lambda1), 1.5) * radius;
}

TEST_F(NavigationTest, NavSetup) {
  Camera cam;
  Quadric quad(-15,0,100,"backprojection");
  Eigen::Vector3d pos = 1e4*Eigen::Vector3d::UnitX();
  Eigen::Vector3d location = Eigen::Vector3d::Zero();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  cam.moveCamera(pos);
  cam.pointTo(location, up_vector);
  Eigen::Matrix4d envelope = quad.getEnvelope();
  std::cout << cam << std::endl;
  Eigen::Matrix3d locus = adjugate(cam.projectQuadric(envelope));
  // Eigen::Matrix3d inv_cam = cam.getInverseIntrinsicMatrix();
  // locus = inv_cam * locus * inv_cam.transpose();
  // Conic conic(locus);
  // std::cout << conic << std::endl;


  Eigen::Matrix3d K = cam.getIntrinsicMatrix();
  Eigen::Matrix3d Kinv = cam.getInverseIntrinsicMatrix();
  Eigen::MatrixXd ext = cam.getExtrinsicMatrix();
  Eigen::MatrixXd prj_mtx = cam.getProjectionMatrix();
  Eigen::Matrix3d proj = prj_mtx * envelope * prj_mtx.transpose();
  ASSERT_TRUE(adjugate(proj).isApprox(locus,1e-12));
  Eigen::Matrix3d symadj = symmetricAdjugate(proj);
  ASSERT_TRUE(symadj.isApprox(locus,1e-12));

  Eigen::Matrix3d exp_envelope = ext * envelope * ext.transpose();
  Eigen::Matrix3d act_envelope = Kinv * proj * Kinv.transpose(); 
  ASSERT_TRUE(exp_envelope.isApprox(act_envelope,1e-12));
  ASSERT_TRUE(K.inverse().isApprox(Kinv));

  Eigen::Matrix3d exp_locus = adjugate(exp_envelope);
  Eigen::Matrix3d act_locus = adjugate(act_envelope); 
  Conic conic(exp_locus);
  std::cout << conic << std::endl;
  std::cout << "Expected:\n" << exp_locus << std::endl;
  // std::cout << "SymAdjugate:\n" << act_locus << std::endl;
  // std::cout << "Actual:\n"   << act_locus << std::endl;
  // std::cout << "Diff:\n"   << act_locus - exp_locus << std::endl;
  ASSERT_TRUE(exp_locus.isApprox(act_locus,1e-8));


  // double smaj = 50, smin = 20, xc = 30, yc = 40, angle = rad2deg(15);
  // Conic conic(smaj, smin, xc, yc, angle);
  // Eigen::Matrix3d locus = conic.getLocus();
  // // Eigen::Matrix3d locus_canonical = canonical(locus);             
  // // // std::cout << "Canonical:\n" << locus_canonical << std::endl;
  // // Conic canon(locus_canonical);
  // // std::cout << canon << std::endl;
  
  Eigen::Vector3d normal;
  double dist;
  conicBackprojection(exp_envelope, quad.getRadius(), normal, dist);
  std::cout << "Quadric radius: " << quad.getRadius() << std::endl;
  std::cout << "Conic normal: " << normal.transpose() << " | distance: " << dist << std::endl;
}

TEST_F(NavigationTest, KanataniEllipseCheck) {
  // double semimajor = 100, semiminor = 40, xc = 40, yc = 120, phi = deg2rad(10);
  // Conic conic(semimajor, semiminor, xc, yc, phi);
}

TEST_F(NavigationTest, ImageConic2PlaneConic) {

}
