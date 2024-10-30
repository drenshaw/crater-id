#include "navigation.h"

#include <eigen3/Eigen/Dense>
#include <vector>
#include <numeric>
#include <iostream>

namespace Shiu {
int findOppositeSignedValueIndex(const std::vector<double>& vec) {
  double sign_val = std::accumulate(
    begin(vec), end(vec), 1.0, std::multiplies<double>());
  auto it = std::find_if(vec.begin(), vec.end(), [sign_val](double num) {
    return (num * sign_val) > 0;
  });
  return it != vec.end() ? it - vec.begin() : -1;
}

double getBackprojectionDistance(const double lambda1, const double lambda2, const double lambda3, const double radius) {
  double abs_l1 = std::abs(lambda1);
  double abs_l2 = std::abs(lambda2);
  double abs_l3 = std::abs(lambda3);
  double scalar = 1/std::abs(lambda2);
  double r = scalar * std::sqrt((abs_l1*abs_l3*(abs_l2 + abs_l3))/(abs_l1 + abs_l3));
  double dist = radius/r;
  return dist;
}

void getBackprojectionLambda1(const std::vector<double>& eigenvalues,
                                  const Eigen::Matrix3d& eigenvectors,
                                  int& mu_d_idx,
                                  Eigen::Vector3d& e3) {
  // Find the eigenvalue with the different sign
  // Eqns 3.1.2.5 - 3.1.2.7
  Eigen::Vector3d f_d;
  mu_d_idx = findOppositeSignedValueIndex(eigenvalues);
  if(mu_d_idx == -1) {
    // std::cerr << "Eigenvalue of opposite sign not found: " << eigenvalues.transpose() << std::endl;
    return;
  }
  // mu_d = eigenvalues.at(mu_d_idx);
  f_d = eigenvectors.col(mu_d_idx);

  // lambda3 = mu_d;
  double dotted = f_d.dot(Eigen::Vector3d::UnitZ());
  assert(dotted != 0);
  e3 = dotted > 0 ? f_d : -f_d;
  // std::cout << "Unique value: " << lambda3 << " at index " << mu_d_idx << std::endl;
  // std::cout << "Remaining eigenvalues: " << eigen_cpy.at(0) << ", " << eigen_cpy.at(1) << std::endl;
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

std::tuple<Eigen::Vector3d, Eigen::Vector3d> getBackprojectionNormals() {

}

void conicBackprojection(const Eigen::Matrix3d& conic, const double radius, Eigen::Vector3d& normal, double& dist) {
  // Applying the concept of the canonical frame from Shiu:
  // "3D loc. of circular and spherical features by monocular ... vision"
  // Eigenvalue decomposition
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(conic);
  if (eigensolver.info() != Eigen::Success) {
    std::cerr << "Eigenvalue decomposition failed!" << std::endl;
    return;
  }
  
  // Eigenvalues and eigenvectors
  Eigen::Vector3d eigenval = eigensolver.eigenvalues();
  Eigen::Matrix3d eigenvec = eigensolver.eigenvectors();
  // Store eigenvalues in a vector for easy manipulation
  std::vector<double> eigenvalues(eigenval.data(), eigenval.data() + eigenval.size());

  double lambda3, lambda2, lambda1;
  Eigen::Vector3d g1, g2, e3, e2, e1;
  int mu_d_idx;  

  // Find the eigenvalue with the different sign
  // Eqns 3.1.2.5 - 3.1.2.7
  getBackprojectionLambda1(eigenvalues, eigenvec, mu_d_idx, e3);
  lambda3 = eigenval(mu_d_idx);
  eigenvalues.erase(eigenvalues.begin() + mu_d_idx);

  // std::cout << "Unique value: " << lambda3 << " at index " << mu_d_idx << std::endl;
  assert(eigenvalues.size() == 2);
  // std::cout << "Remaining eigenvalues: " << eigenvalues.at(0) << ", " << eigenvalues.at(1) << std::endl;
  assert(eigenvalues.at(0) * eigenvalues.at(1) > 0);
  
  // Get remaining lambdas
  // Eqns 3.1.2.8 - 3.1.2.11
  getBackprojectionRemainingLambdas(eigenvalues, eigenvec, e3, lambda1, lambda2, e1, e2);

  // Eqn 3.1.3.6
  dist = getBackprojectionDistance(lambda1, lambda2, lambda3, radius);
}

} // end namespace
