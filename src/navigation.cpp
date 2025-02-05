#include "navigation.h"
#include "quadrics.h"
#include "math_utils.h"

#include <eigen3/Eigen/Dense>
#include <vector>
#include <numeric>
#include <iostream>

/*********************************************************/
/***************  Single Crater Functions  ***************/
/*********************************************************/
namespace Preprocessing {

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

void getEigenParts( const Eigen::Matrix3d& conic, 
                    double& lambda1, double& lambda2, double& lambda3, 
                    Eigen::Vector3d& u1, Eigen::Vector3d& u2, Eigen::Vector3d& u3) {

  // Eigenvalues and eigenvectors
  Eigen::Vector3d eigenval;
  Eigen::Matrix3d eigenvec;
  if(!getEigenstuffConic(conic, eigenval, eigenvec)) {
    std::cerr << __func__ << "-> backprojection eigenstuff failed.\n";
    return;
  }
  int mu_d_idx;  

  // Find the eigenvalue with the different sign
  // Eqns 3.1.2.5 - 3.1.2.7
  std::vector<int> indices = {0, 1, 2};
  Preprocessing::getBackprojectionLambda3(eigenval, eigenvec, mu_d_idx, u3);
  lambda3 = eigenval(mu_d_idx);
  // std::cout << "Eigenvalues: " << eigenval.transpose() << std::endl;
  // std::cout << "Eigenvectors:\n" << eigenvec << std::endl;
  indices.erase(indices.begin() + mu_d_idx);
  Eigen::Vector2d rem_eigenval = eigenval(indices);
  Eigen::MatrixXd rem_eigenvec = eigenvec(Eigen::all, indices);

  // Get remaining lambdas
  // Eqns 3.1.2.8 - 3.1.2.11
  Preprocessing::getBackprojectionLambda2(rem_eigenval, rem_eigenvec, lambda1, lambda2, u2);
  u1 = u2.cross(u3);
}

void getBackprojectionLambda3(const Eigen::Vector3d& eigenvalues,
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

} // namespace Preprocessing

namespace Christian {

void getBackprojectedCenter(const double radius, const std::array<double, 3>& lambdas,
                            const Eigen::Matrix3d& canonizing, 
                            std::array<Eigen::Vector3d, 2>& centers) {

  double lambda1 = lambdas.at(0);
  double lambda2 = lambdas.at(1);
  double lambda3 = lambdas.at(2);
  double kx2 = 1/std::abs(lambda1);
  double ky2 = 1/std::abs(lambda2);
  double kz2 = 1/std::abs(lambda3);
  double alpha = (ky2 - kx2) / (ky2 + kz2);
  double kxz2 = kx2/kz2;
  double sqrtAlphaPlusKxz = std::sqrt(alpha + kxz2);
  Eigen::Vector3d center_posX, center_negX, center_pos, center_neg;
  double scalarR = radius/sqrtAlphaPlusKxz;
  center_posX <<  std::sqrt(alpha*kxz2), 0, 1;
  center_negX << -std::sqrt(alpha*kxz2), 0, 1;
  center_posX *= scalarR;
  center_negX *= scalarR;
  center_pos  = canonizing * center_posX;
  center_neg  = canonizing * center_negX;
  centers = {center_pos, center_neg};
}

void getBackprojectedNormal(const std::array<double, 3>& lambdas,
                            const Eigen::Matrix3d& canonizing,
                            std::array<Eigen::Vector3d, 2>& normals) {

  double lambda1 = lambdas.at(0);
  double lambda2 = lambdas.at(1);
  double lambda3 = lambdas.at(2);  
  // X-axis corresponds to semiminor (smaller) axis
  double kx2 = 1/std::abs(lambda1);
  double ky2 = 1/std::abs(lambda2);
  double kz2 = 1/std::abs(lambda3);
  double alpha = (ky2 - kx2) / (ky2 + kz2);        
  Eigen::Vector3d normal_posX, normal_negX, normal_pos, normal_neg;
  double kxz2 = kx2/kz2;
  double kxz = std::sqrt(kxz2);
  double sqrtAlphaPlusKxz = std::sqrt(alpha + kxz2);
  double scalarN = 1/sqrtAlphaPlusKxz;
  normal_posX <<  std::sqrt(alpha), 0, -kxz;
  normal_negX << -std::sqrt(alpha), 0, -kxz;
  normal_posX *= scalarN;
  normal_negX *= scalarN;
  normal_pos  = canonizing * normal_posX;
  normal_neg  = canonizing * normal_negX;
  normals = {normal_pos, normal_neg};
}

void conicBackprojection( const Eigen::Matrix3d& conic_locus, const double radius, 
                          std::array<Eigen::Vector3d, 2>& centers, 
                          std::array<Eigen::Vector3d, 2>& normals) {
  // Applying the transforms found in Christian white paper
  // "Perspective projection of ellipses and ellipsoids with applications to spacecraft navigation"
  double lambda1, lambda2, lambda3;
  Eigen::Vector3d u1, u2, u3;
  Preprocessing::getEigenParts(conic_locus, lambda1, lambda2, lambda3, u1, u2, u3);
  Eigen::Matrix3d canonizing;
  canonizing << u1, u2, u3;

  std::array<double, 3> lambdas = {lambda1, lambda2, lambda3};
  getBackprojectedCenter(radius, lambdas, canonizing, centers);
  getBackprojectedNormal(lambdas, canonizing, normals);
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

void conicBackprojection( const Eigen::Matrix3d& conic_locus, const double radius, 
                              std::array<Eigen::Vector3d, 2>& centers, 
                              std::array<Eigen::Vector3d, 2>& normals) {
  // Applying the concept of the canonical frame from Shiu:
  // "3D loc. of circular and spherical features by monocular ... vision"
  // Eigenvalues and eigenvectors
  double lambda1, lambda2, lambda3;
  Eigen::Vector3d u1, u2, u3;
  Preprocessing::getEigenParts(conic_locus, lambda1, lambda2, lambda3, u1, u2, u3);
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

double getDistanceFromConicLocus(const Eigen::Matrix3d& conic_locus, const double radius) {
  double lambda1, lambda2, lambda3;
  Eigen::Vector3d u1, u2, u3;
  Preprocessing::getEigenParts(conic_locus, lambda1, lambda2, lambda3, u1, u2, u3);
  double dist = getBackprojectionDistance(lambda1, lambda2, lambda3, radius);
  return dist;
}

} // end namespace Shiu

namespace Kanatani {
  // TODO: the paper referenced ("3D interpretation of conics and orthogonality"), only
  // uses one of two possible orientations and centers. Figure out how to use the same
  // framework to get the other orientation and center.
void conicBackprojection( const Eigen::Matrix3d& conic_locus, const double radius, 
                          Eigen::Vector3d& normal, double& dist) {
  // Applying the concept of the canonical frame from Kanatani:
  // "3D interpretation of conics and orthogonality"
  double lambda1, lambda2, lambda3;
  Eigen::Vector3d u1, u2, u3;
  // We need to switch the 1 & 2 indices for Kanatani's method
  Preprocessing::getEigenParts(conic_locus, lambda2, lambda1, lambda3, u2, u1, u3);

  double scalar1 = std::sqrt((lambda2 - lambda1)/(lambda2 - lambda3));
  double scalar2 = std::sqrt((lambda1 - lambda3)/(lambda2 - lambda3));
  assert(!std::isnan(scalar1));
  assert(!std::isnan(scalar2));
  // TODO: do I need this negative sign?
  normal = -(scalar1*u2 + scalar2*u3);
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

/*********************************************************/
/*************** Multiple Crater Functions ***************/
/*********************************************************/
uint chooseSupportingPlanes(const double angle, 
                            const std::array<Eigen::Vector3d, 2>& normals1, 
                            const std::array<Eigen::Vector3d, 2>& normals2) {
  // Kanatani calls the plane in the camera frame containing the 3d circle
  // the "supporting plane"
  double angle00 = getAngleBetweenVectors(normals1.at(0), normals2.at(0));
  double angle01 = getAngleBetweenVectors(normals1.at(0), normals2.at(1));
  double angle10 = getAngleBetweenVectors(normals1.at(1), normals2.at(0));
  double angle11 = getAngleBetweenVectors(normals1.at(1), normals2.at(1));
  // TODO: Getting issues when the cone is circular (i.e., crater is perpendicular to us)
  Eigen::Vector4d angles;
  angles << angle00, angle01, angle10, angle11;
  Eigen::Vector4d d_angles = angles.array() - angle;
  Eigen::Vector4d absd_angles = d_angles.cwiseAbs();
  Eigen::Vector4d::Index min_idx;
  absd_angles.minCoeff(&min_idx);
  const double max_angle = rad2deg(3.);
  if(absd_angles(min_idx) > max_angle) {
    // std::cout << __func__ << "-> Truth angle: " << rad2deg(angle)
    //           << " vs " << rad2deg(angles(min_idx))
    //           << " @ " << min_idx  << std::endl;
    // std::cout << "\tAll angles: " << angles.transpose() << std::endl;
    std::cout << "Angle exceeds max angle" << max_angle 
              << " : " << d_angles(min_idx) 
              << " @ index " << min_idx << " | " 
              << d_angles.transpose() << std::endl;
    // TODO: This is a great choice for std::optional
    return 5;
  }
  return min_idx;
}

void selectSupportingPlaneIndex(const uint index, uint& index_a, uint& index_b) {
  index_a = index / 2;
  index_b = index % 2;
}

void selectSupportingPlaneNormals(const uint index,
                                  const std::array<Eigen::Vector3d, 2>& normals1,
                                  const std::array<Eigen::Vector3d, 2>& normals2,
                                  Eigen::Vector3d& normal1,
                                  Eigen::Vector3d& normal2) {
  uint index_a, index_b;
  selectSupportingPlaneIndex(index, index_a, index_b);
  normal1 = normals1.at(index_a);
  normal2 = normals2.at(index_b);
}

void selectSupportingPlaneCenters(const uint index,
                                  const std::array<Eigen::Vector3d, 2>& centers1,
                                  const std::array<Eigen::Vector3d, 2>& centers2,
                                  Eigen::Vector3d& center1,
                                  Eigen::Vector3d& center2) {
  selectSupportingPlaneNormals(index, centers1, centers2, center1, center2);
}

void reprojectLociiToQuadrics(const std::vector<Quadric>& quadrics,
                              const std::vector<Eigen::Matrix3d>& locii,
                              std::vector<Eigen::Vector3d>& centers,
                              std::vector<Eigen::Vector3d>& normals) {
  uint q_size = quadrics.size();
  centers.clear();
  normals.clear();
  std::vector<uint> indices;
  indices.reserve(q_size - 1);
  std::vector<Eigen::Matrix3d>::const_iterator it;
  for (it = locii.begin(); it != locii.end(); it++) {
    int index = getIndex(locii.begin(), it);
    Quadric quad = quadrics.at(index);
    double radius1 = quad.getRadius();
    std::array<Eigen::Vector3d, 2> centers1, normals1;
    Christian::conicBackprojection(*it, radius1, centers1, normals1);

    // Now remove the current quadric and locus; put them in separate vector
    std::vector<Quadric> other_quadrics = copyAllBut(quadrics, quad);
    std::vector<Eigen::Matrix3d> other_locii = copyAllBut(locii, *it);
    Eigen::Vector3d center1, center2, normal1, normal2;

    std::vector<Eigen::Matrix3d>::iterator it_other;
    const double max_angle = deg2rad(3.);
    for(it_other = other_locii.begin(); it_other != other_locii.end(); it_other++) {
      int index_other = getIndex(other_locii.begin(), it_other);
      Quadric other_quadric = other_quadrics.at(index_other);
      double radius2 = other_quadric.getRadius();
      std::array<Eigen::Vector3d, 2> centers2, normals2;
      Christian::conicBackprojection(*it_other, radius2, centers2, normals2);
      double angle = quad.getAngleBetweenQuadrics(other_quadric);
      uint idx_ab = chooseSupportingPlanes(angle, normals1, normals2);
      selectSupportingPlaneNormals(idx_ab, normals1, normals2, normal1, normal2);
      selectSupportingPlaneCenters(idx_ab, centers1, centers2, center1, center2);
      double angle_calc = getAngleBetweenVectors(normal1, normal2);
      if(std::abs(angle_calc - angle)>max_angle) {
        std::cerr << "Angle " << angle << " exceeds max angle: " << max_angle << std::endl;
      }
      // assert(std::abs(angle_calc - angle)<max_angle);
      indices.push_back(idx_ab / 2);
    }
    if(!allEqualVector(indices)) {
      std::cout << "Not all entries are equal.\n";
    }
    normals.push_back(normal1);
    centers.push_back(center1);
    indices.clear();
  }
}

void reprojectionsToPlanes( const std::vector<Eigen::Vector3d>& centers,
                            const std::vector<Eigen::Vector3d>& normals,
                            std::vector<Eigen::Hyperplane<double, 3> >& planes) {
  std::vector<Eigen::Vector3d>::const_iterator the_chosen;
  for(the_chosen = normals.begin(); the_chosen != normals.end(); the_chosen++) {
    int index = getIndex(normals.begin(), the_chosen);
    Eigen::Vector3d center = centers.at(index);
    Eigen::Hyperplane<double, 3> plane(*the_chosen, center);
    planes.push_back(plane);
  }
}

void calculateHomography( const std::vector<Eigen::Hyperplane<double, 3> >& planes_world,
                          const std::vector<Eigen::Hyperplane<double, 3> >& planes_cam,
                          Eigen::Quaterniond& attitude, Eigen::Vector3d& position) {
  const int min_n_craters = 3;
  if(planes_world.size() < min_n_craters) {
    std::cerr << "Not enough craters to perform estimation.\n";
    return;
  }                            
  uint PLANE_DIM = 4;
  uint n_planes = planes_world.size();
  assert(n_planes == planes_cam.size());
  Eigen::VectorXd lhs(PLANE_DIM*n_planes);
  Eigen::MatrixXd rhs(PLANE_DIM*n_planes, PLANE_DIM*(PLANE_DIM-1));
  for(size_t i = 0; i < n_planes; i++) {
    lhs.block<4,1>(4*i,0) = planes_cam.at(i).coeffs();
    Eigen::Hyperplane<double, 3> world_plane = planes_world.at(i);
    lhs(PLANE_DIM*i+PLANE_DIM-1) -= world_plane.offset();
    
    double a = world_plane.normal()(0);
    Eigen::Matrix4d diag_a = a*Eigen::Matrix4d::Identity();
    double b = world_plane.normal()(1);
    Eigen::Matrix4d diag_b = b*Eigen::Matrix4d::Identity();
    double c = world_plane.normal()(2);
    Eigen::Matrix4d diag_c = c*Eigen::Matrix4d::Identity();
    rhs.block<4,4>(PLANE_DIM*i, PLANE_DIM*0) = diag_a;
    rhs.block<4,4>(PLANE_DIM*i, PLANE_DIM*1) = diag_b;
    rhs.block<4,4>(PLANE_DIM*i, PLANE_DIM*2) = diag_c;
  }

  Eigen::BDCSVD<Eigen::MatrixXd> SVD(rhs, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd solu = SVD.solve(lhs);
  Eigen::MatrixXd P_est = solu.reshaped(4,3).transpose();
  Eigen::Matrix3d dcm = P_est.topLeftCorner(3,3);
  attitude = Eigen::Quaterniond(dcm).inverse();
  position = P_est.topRightCorner(3,1);
}

void solve_navigation_problem(const std::vector<Quadric>& quadrics,
                              const std::vector<Eigen::Matrix3d>& locii,
                              Eigen::Quaterniond& attitude, Eigen::Vector3d& position) {
  std::vector<Eigen::Hyperplane<double, 3> > planes_world;
  for(std::vector<Quadric>::const_iterator it = quadrics.begin(); it != quadrics.end(); it++) {
    planes_world.push_back((*it).getPlane());
  }
  std::vector<Eigen::Vector3d> centers;
  std::vector<Eigen::Vector3d> normals;
  reprojectLociiToQuadrics(quadrics, locii, centers, normals);
  std::vector<Eigen::Hyperplane<double, 3> > planes_cam;
  reprojectionsToPlanes(centers, normals, planes_cam);
  calculateHomography(planes_world, planes_cam, attitude, position);
}

double calcAngleFromEllipseAxes(const double Rmax, const double Rmin, const double f) {
  // Using Eq 16 from "Constraints on quadratic curves under perspective projection"
  // by Safaee-Rad, Smith, Benhabib, and Tchoukanov
  // NOTE: the corrected equation uses a minus sign in the denominator
  assert(Rmax >= Rmin && Rmin >= 0);
  const double f_Rmin_2 = std::pow(f/Rmin, 2);
  const double f_Rmax_2 = std::pow(f/Rmax, 2);
  const double num = f_Rmin_2 - f_Rmax_2;
  const double den = f_Rmin_2 - 1;
  const double sin_angle = std::sqrt(num/den);
  return std::asin(sin_angle);
}

double calcAngleFromEllipseAxes(const Eigen::Matrix3d& conic_locus) {
  // Using Eq 3.1.3.9 from "3D location of circular and spherical features by
  // monocular model-based vision" by Shiu and Ahmed
  double lambda1, lambda2, lambda3;
  Eigen::Vector3d u1, u2, u3;
  Preprocessing::getEigenParts(conic_locus, lambda1, lambda2, lambda3, u1, u2, u3);
  const double abs_l1 = std::abs(lambda1);
  const double abs_l2 = std::abs(lambda2);
  const double abs_l3 = std::abs(lambda3);
  const double radical = (abs_l1 - abs_l2)/(abs_l2 + abs_l3);
  const double phi = std::atan(std::sqrt(radical));
  return phi;
}

double attitudeError(const Eigen::Quaterniond& Qest, const Eigen::Quaterniond& Qtrue) {
  Eigen::Quaterniond Qdiff = Qest.normalized().inverse() * Qtrue.normalized();
  Qdiff.normalize();
  return wrap_2pi(2*std::acos(Qdiff.w()));
}

