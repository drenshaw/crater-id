#include "quadrics.h"

#include <eigen3/Eigen/Geometry>
#include <numeric>


/*********************************************************/
/***************  Single Crater Functions  ***************/
/*********************************************************/

namespace Preprocessing {
int findOppositeSignedValueIndex(const std::vector<double>& vec);
bool getEigenstuffConic(const Eigen::Matrix3d& conic, Eigen::Vector3d& eigenval, Eigen::Matrix3d& eigenvec);

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
void getEigenParts( const Eigen::Matrix3d& conic, 
                    double& lambda1, double& lambda2, double& lambda3, 
                    Eigen::Vector3d& u1, Eigen::Vector3d& u2, Eigen::Vector3d& u3);
double getBackprojectionDistance( const double lambda1, 
                                  const double lambda2, 
                                  const double lambda3, 
                                  const double radius);
void getBackprojectionLambda3(const Eigen::Vector3d& eigenvalues,
                              const Eigen::Matrix3d& eigenvectors,
                              int& mu_d_idx,
                              Eigen::Vector3d& e3);
void getBackprojectionLambda2(const Eigen::Vector2d& eigenvalues,
                              const Eigen::MatrixXd& eigenvectors,
                              double& lambda1, double& lambda2,
                              Eigen::Vector3d& e2);
}

namespace Christian {
                               
void conicBackprojection( const Eigen::Matrix3d& conic_locus, const double radius, 
                              std::array<Eigen::Vector3d, 2>& centers, 
                              std::array<Eigen::Vector3d, 2>& normals);
void getBackprojectedCenter(const double radius, const std::array<double, 3>& lambdas,
                            const Eigen::Matrix3d& canonizing, 
                            std::array<Eigen::Vector3d, 2>& centers);

void getBackprojectedNormal(const std::array<double, 3>& lambdas,
                            const Eigen::Matrix3d& canonizing,
                            std::array<Eigen::Vector3d, 2>& normals);
} // end namespace Christian

namespace Shiu {
void getBackprojectionNormalCanonical(const double lambda1, const double lambda2, 
                                          const double lambda3, Eigen::Vector3d& normal);
void getBackprojectionNormal( const double lambda1, const double lambda2, const double lambda3,
                                  const Eigen::Vector3d& e1, Eigen::Vector3d e3, 
                                  std::array<Eigen::Vector3d, 2>& normals);
void getBackprojectionCenterCanonical(
  const double radius, const double lambda1, const double lambda2, const double lambda3, Eigen::Vector3d& center_offset);
void getBackprojectionCenter( const double radius, const double lambda1, const double lambda2, const double lambda3,
                                  const Eigen::Vector3d& e1, Eigen::Vector3d e3, 
                                  std::array<Eigen::Vector3d, 2> centers);
void conicBackprojection( const Eigen::Matrix3d& conic_locus, const double radius, 
                              std::array<Eigen::Vector3d, 2>& centers, 
                              std::array<Eigen::Vector3d, 2>& normals);
double getDistanceFromConicLocus(const Eigen::Matrix3d& conic_locus, const double radius);
} // end namespace Shiu

namespace Kanatani {
void conicBackprojection( const Eigen::Matrix3d& conic_locus, const double radius, 
                          Eigen::Vector3d& normal, double& dist);
Eigen::Matrix3d canonical(const Eigen::Matrix3d& image_conic);
}

/*********************************************************/
/*************** Multiple Crater Functions ***************/
/*********************************************************/
uint chooseSupportingPlanes(const double angle, 
                            const std::array<Eigen::Vector3d, 2>& normals1, 
                            const std::array<Eigen::Vector3d, 2>& normals2);

void selectSupportingPlaneIndex(const uint index, uint& index_a, uint& index_b);

void selectSupportingPlaneNormals(const uint index,
                                  const std::array<Eigen::Vector3d, 2>& normals1,
                                  const std::array<Eigen::Vector3d, 2>& normals2,
                                  Eigen::Vector3d& normal1,
                                  Eigen::Vector3d& normal2);

void selectSupportingPlaneCenters(const uint index,
                                  const std::array<Eigen::Vector3d, 2>& centers1,
                                  const std::array<Eigen::Vector3d, 2>& centers2,
                                  Eigen::Vector3d& center1,
                                  Eigen::Vector3d& center2);

void reprojectLociiToQuadrics(const std::vector<Quadric>& quadrics,
                              const std::vector<Eigen::Matrix3d>& locii,
                              std::vector<Eigen::Vector3d>& centers,
                              std::vector<Eigen::Vector3d>& normals);


void reprojectionsToPlanes( const std::vector<Eigen::Vector3d>& centers,
                            const std::vector<Eigen::Vector3d>& normals,
                            std::vector<Eigen::Hyperplane<double, 3> >& planes);

void calculateHomography( const std::vector<Eigen::Hyperplane<double, 3> >& planes_world,
                          const std::vector<Eigen::Hyperplane<double, 3> >& planes_cam,
                          Eigen::Quaterniond& attitude, Eigen::Vector3d& position);

void solve_navigation_problem(const std::vector<Quadric>& quadrics,
                              const std::vector<Eigen::Matrix3d>& locii,
                              Eigen::Quaterniond& attitude, Eigen::Vector3d& position);
double calcAngleFromEllipseAxes(const double Rmax, const double Rmin, const double f);
double calcAngleFromEllipseAxes(const Eigen::Matrix3d& conic_locus);
double attitudeError(const Eigen::Quaterniond& Qest, const Eigen::Quaterniond& Qtrue);