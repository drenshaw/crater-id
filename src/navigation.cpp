#include "navigation.h"

#include <eigen3/Eigen/Dense>
#include <vector>
#include <numeric>
#include <iostream>

int findOppositeSignedValueIndex(const std::vector<double>& vec) {
  double sign_val = std::accumulate(
    begin(vec), end(vec), 1.0, std::multiplies<double>());
  auto it = std::find_if(vec.begin(), vec.end(), [sign_val](double num) {
    return (num * sign_val) > 0;
  });
  return it != vec.end() ? it - vec.begin() : -1;
}

bool getEigenstuffConic(const Eigen::Matrix3d& conic, Eigen::Vector3d& eigenval, Eigen::Matrix3d& eigenvec) {
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(conic);
  if (eigensolver.info() != Eigen::Success) {
    std::cerr << "Eigenvalue decomposition failed!" << std::endl;
    return false;
  }
  
  // Eigenvalues and eigenvectors
  eigenval = eigensolver.eigenvalues();
  eigenvec = eigensolver.eigenvectors();
  return true;
}

namespace Christian {

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
  Eigen::Vector2d::Index maxRow,maxCol;
  eig.maxCoeff(&maxRow, &maxCol);
  std::cout << "Max row: " << maxRow << ", max col: " << maxCol << std::endl;
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
  std::vector<int> indices = {0, 1, 2};
  Shiu::getBackprojectionLambda1(eigenval, eigenvec, mu_d_idx, u3);
  lambda3 = eigenval(mu_d_idx);
  indices.erase(indices.begin() + mu_d_idx);
  Eigen::Vector2d rem_eigenval = eigenval(indices);
  Eigen::MatrixXd rem_eigenvec = eigenvec(Eigen::all, indices);

  // Get remaining lambdas
  // Eqns 3.1.2.8 - 3.1.2.11
  Shiu::getBackprojectionLambda2(rem_eigenval, rem_eigenvec, lambda1, lambda2, u2);
  u1 = u2.cross(u3);
  // Eqn 3.1.3.6
  // double dist = getBackprojectionDistance(lambda1, lambda2, lambda3, radius);
  Shiu::getBackprojectionCenter(radius, lambda1, lambda2, lambda3, u1, u3, centers);
  Shiu::getBackprojectionNormal(lambda1, lambda2, lambda3, u1, u3, normals);
  Eigen::Vector3d center1, center2, normal1, normal2;
  center1 = centers.at(0);
  center2 = centers.at(1);
  normal1 = normals.at(0);
  normal2 = normals.at(1);
  std::cout << "Center1: " << center1.transpose() << " | Distance: " << center1.norm() << std::endl;
  std::cout << "Normal1: " << normal1.transpose() << std::endl;
  std::cout << "Center2: " << center2.transpose() << " | Distance: " << center2.norm() << std::endl;
  std::cout << "Normal2: " << normal2.transpose() << std::endl;
}

} // end namespace Christian

namespace Shiu {

double getBackprojectionDistance(const double lambda1, const double lambda2, const double lambda3, const double radius) {
  // Calculates the distance from the camera to the backprojected ellipse center
  // Note: this is generally not the true quadric center
  double abs_l1 = std::abs(lambda1);
  double abs_l2 = std::abs(lambda2);
  double abs_l3 = std::abs(lambda3);
  double scalar = 1/std::abs(lambda2);
  double r = scalar * std::sqrt((abs_l1*abs_l3*(abs_l2 + abs_l3))/(abs_l1 + abs_l3));
  double dist = radius/r;
  return dist;
}

void getBackprojectionLambda1(const Eigen::Vector3d& eigenvalues,
                              const Eigen::Matrix3d& eigenvectors,
                              int& mu_d_idx,
                              Eigen::Vector3d& e3) {
  // Find the eigenvalue with the different sign
  // Shiu Eqns 3.1.2.5 - 3.1.2.7
  // Christian Eqn 58
  Eigen::Vector3d f_d;
  std::vector<double> eigenval(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
  mu_d_idx = findOppositeSignedValueIndex(eigenval);
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

void getBackprojectionLambda2(const Eigen::Vector2d& eigenvalues,
                              const Eigen::MatrixXd& eigenvectors,
                              double& lambda1, double& lambda2,
                              Eigen::Vector3d& e2) {
  assert(eigenvectors.cols() == 2);
  Eigen::Vector3d g1, g2;
  // Get remaining lambdas
  // Eqns 3.1.2.8 - 3.1.2.11
  
  Eigen::Vector2d::Index minRow, maxRow;
  Eigen::Vector2d abs_eigenvalues = eigenvalues.cwiseAbs();
  abs_eigenvalues.minCoeff(&minRow);
  abs_eigenvalues.maxCoeff(&maxRow);
  lambda1 = eigenvalues(maxRow);
  lambda2 = eigenvalues(minRow);
  e2 = eigenvectors.col(minRow);
  // double omega1, omega2;
  // // e1 = eigenvectors.col(minRow);
  // if(omega1 == omega2) {
  //   std::cerr << "Not confident in handling right cones yet!\n";
  // }
  // if(std::abs(omega1) < std::abs(omega2)) {
  //   lambda2 = omega1;
  //   lambda1 = omega2;
  //   e2 = g1;
  // }
  // else {
  //   lambda2 = omega2;
  //   lambda1 = omega1;
  //   e2 = g2;
  // }
}

void getBackprojectionNormalCanonical(const double lambda1, const double lambda2, 
                                      const double lambda3, Eigen::Vector3d& normal) {
  double abs_l1 = std::abs(lambda1);
  double abs_l2 = std::abs(lambda2);
  double abs_l3 = std::abs(lambda3);
  double n_x =  std::sqrt((abs_l1-abs_l2)/(abs_l1+abs_l3));
  double n_y =  0;
  double n_z = -std::sqrt((abs_l2+abs_l3)/(abs_l1+abs_l3));
  normal << n_x, n_y, n_z;
}

void getBackprojectionNormal( const double lambda1, const double lambda2, const double lambda3,
                              const Eigen::Vector3d& e1, Eigen::Vector3d e3, 
                              std::array<Eigen::Vector3d, 2>& normals) {
  Eigen::Vector3d canonical_normal;
  getBackprojectionNormalCanonical(lambda1, lambda2, lambda3, canonical_normal);
  double e1x = e1(0);
  double e1y = e1(1);
  double e1z = e1(2);
  double e3x = e3(0);
  double e3y = e3(1);
  double e3z = e3(2);
  double n_x = canonical_normal(0);
  double n_z = canonical_normal(2);
  double n_x1 =  e1x*n_x + e3x*n_z;
  double n_y1 =  e1y*n_x + e3y*n_z;
  double n_z1 =  e1z*n_x + e3z*n_z;
  double n_x2 = -e1x*n_x + e3x*n_z;
  double n_y2 = -e1y*n_x + e3y*n_z;
  double n_z2 = -e1z*n_x + e3z*n_z;
  Eigen::Vector3d normal1, normal2;
  normal1 << n_x1, n_y1, n_z1;
  normal2 << n_x2, n_y2, n_z2;
  normals = {normal1, normal2};
}

void getBackprojectionCenterCanonical(
  const double radius, const double lambda1, const double lambda2, const double lambda3, Eigen::Vector3d& center_offset) {
  // Eqn 3.1.3.12 - 3.1.3.15
  double abs_l1 = std::abs(lambda1);
  double abs_l2 = std::abs(lambda2);
  double abs_l3 = std::abs(lambda3);
  double x = radius * std::sqrt((abs_l3*(abs_l1-abs_l2))/(abs_l1*(abs_l1+abs_l3)));
  double y = 0;
  double z = radius * std::sqrt((abs_l1*(abs_l2+abs_l3))/(abs_l3*(abs_l1+abs_l3)));
  center_offset << x, y, z;
}

void getBackprojectionCenter( const double radius, const double lambda1, const double lambda2, const double lambda3,
                              const Eigen::Vector3d& e1, Eigen::Vector3d e3, 
                              std::array<Eigen::Vector3d, 2> centers) {
  Eigen::Vector3d canonical_center;
  getBackprojectionCenterCanonical(radius, lambda1, lambda2, lambda3, canonical_center);
  double e1x = e1(0);
  double e1y = e1(1);
  double e1z = e1(2);
  double e3x = e3(0);
  double e3y = e3(1);
  double e3z = e3(2);
  double n_x = canonical_center(0);
  double n_z = canonical_center(2);
  double n_x1 =  e1x*n_x + e3x*n_z;
  double n_y1 =  e1y*n_x + e3y*n_z;
  double n_z1 =  e1z*n_x + e3z*n_z;
  double n_x2 = -e1x*n_x + e3x*n_z;
  double n_y2 = -e1y*n_x + e3y*n_z;
  double n_z2 = -e1z*n_x + e3z*n_z;
  Eigen::Vector3d center1, center2;
  center1 << n_x1, n_y1, n_z1;
  center2 << n_x2, n_y2, n_z2;
  centers = {center1, center2};
}                                    

void conicBackprojection( const Eigen::Matrix3d& conic, const double radius, 
                              std::array<Eigen::Vector3d, 2>& centers, 
                              std::array<Eigen::Vector3d, 2>& normals) {
  // Applying the concept of the canonical frame from Shiu:
  // "3D loc. of circular and spherical features by monocular ... vision"
  // Eigenvalues and eigenvectors
  Eigen::Vector3d eigenvalues;
  Eigen::Matrix3d eigenvectors;
  if(!getEigenstuffConic(conic, eigenvalues, eigenvectors)) {
    std::cerr << __func__ << "-> backprojection eigenstuff failed.\n";
    return;
  }

  double lambda3, lambda2, lambda1;
  Eigen::Vector3d g1, g2, u3, u2, u1;
  int mu_d_idx;  

  // Find the eigenvalue with the different sign (--+ or -++)
  // Eqns 3.1.2.5 - 3.1.2.7
  std::vector<int> indices = {0, 1, 2};
  getBackprojectionLambda1(eigenvalues, eigenvectors, mu_d_idx, u3);
  lambda3 = eigenvalues(mu_d_idx);
  std::cout << "Eigenvalues: " << eigenvalues.transpose() << std::endl;
  std::cout << "Eigenvectors:\n" << eigenvectors << std::endl;
  indices.erase(indices.begin() + mu_d_idx);
  Eigen::Vector2d rem_eigenval = eigenvalues(indices);
  Eigen::MatrixXd rem_eigenvec = eigenvectors(Eigen::all, indices);

  // Get remaining lambdas
  // Eqns 3.1.2.8 - 3.1.2.11
  getBackprojectionLambda2(rem_eigenval, rem_eigenvec, lambda1, lambda2, u2);
  u1 = u2.cross(u3);
  // Eqn 3.1.3.6
  // double dist = getBackprojectionDistance(lambda1, lambda2, lambda3, radius);
  Shiu::getBackprojectionCenter(radius, lambda1, lambda2, lambda3, u1, u3, centers);
  Shiu::getBackprojectionNormal(lambda1, lambda2, lambda3, u1, u3, normals);
  Eigen::Vector3d center1, center2, normal1, normal2;
  center1 = centers.at(0);
  center2 = centers.at(1);
  normal1 = normals.at(0);
  normal2 = normals.at(1);
  std::cout << "Center1: " << center1.transpose() << " | Distance: " << center1.norm() << std::endl;
  std::cout << "Normal1: " << normal1.transpose() << std::endl;
  std::cout << "Center2: " << center2.transpose() << " | Distance: " << center2.norm() << std::endl;
  std::cout << "Normal2: " << normal2.transpose() << std::endl;
}
} // end namespace Shiu

namespace Kanatani {
void conicBackprojection( const Eigen::Matrix3d& conic_locus, const double radius, 
                          Eigen::Vector3d& normal, double& dist) {
  // Applying the concept of the canonical frame from Kanatani:
  // "3D interpretation of conics and orthogonality"
  // Eigenvalues and eigenvectors
  Eigen::Vector3d eigenvalues;
  Eigen::Matrix3d eigenvectors;
  if(!getEigenstuffConic(conic_locus, eigenvalues, eigenvectors)) {
    std::cerr << __func__ << "-> backprojection eigenstuff failed.\n";
    return;
  }
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
  // assert(lambda3 < 0);
  // assert(lambda1 > 0);
  // assert(lambda3 < lambda2);
  // assert(lambda1 < lambda2);
  double scalar1 = std::sqrt((lambda2 - lambda1)/(lambda2 - lambda3));
  double scalar2 = std::sqrt((lambda1 - lambda3)/(lambda2 - lambda3));
  assert(!std::isnan(scalar1));
  assert(!std::isnan(scalar2));
  normal = scalar1*u2 + scalar2*u3;
  std::cout << "Quadric radius: " << radius << std::endl;
  dist = std::pow(std::abs(lambda1), 1.5) * radius;
}

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

} // end namespace Kanatani

