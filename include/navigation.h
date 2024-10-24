#include <eigen3/Eigen/Geometry>


namespace Shiu {
int findOppositeSignedValueIndex(const std::vector<double>& vec);
double getBackprojectionDistance(const double lambda1, const double lambda2, const double lambda3, const double radius);
void getBackprojectionLambda1(const std::vector<double>& eigenvalues,
                                  const Eigen::Matrix3d& eigenvectors,
                                  int& mu_d_idx,
                                  Eigen::Vector3d& e3);
void getBackprojectionRemainingLambdas( const std::vector<double>& eigenvalues,
                                            const Eigen::Matrix3d& eigenvectors,
                                            const Eigen::Vector3d e3,
                                            double& lambda1, double& lambda2,
                                            Eigen::Vector3d& e1, Eigen::Vector3d& e2);
void conicBackprojection(const Eigen::Matrix3d& conic, const double radius, Eigen::Vector3d& normal, double& dist);
std::tuple<Eigen::Vector3d, Eigen::Vector3d> getBackprojectionNormals();
} // end namespace
