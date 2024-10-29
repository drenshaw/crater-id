#include <eigen3/Eigen/Geometry>
#include <numeric>


int findOppositeSignedValueIndex(const std::vector<double>& vec);

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
bool getEigenstuffConic(const Eigen::Matrix3d& conic, Eigen::Vector3d& eigenval, Eigen::Matrix3d& eigenvec);

namespace Christian {
                               
void conicBackprojection( const Eigen::Matrix3d& conic, const double radius, 
                              std::array<Eigen::Vector3d, 2>& centers, 
                              std::array<Eigen::Vector3d, 2>& normals);
}
namespace Shiu {
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
void conicBackprojection( const Eigen::Matrix3d& conic, const double radius, 
                              std::array<Eigen::Vector3d, 2>& centers, 
                              std::array<Eigen::Vector3d, 2>& normals);
} // end namespace Shiu

namespace Kanatani {
void conicBackprojection( const Eigen::Matrix3d& conic_locus, const double radius, 
                          Eigen::Vector3d& normal, double& dist);
Eigen::Matrix3d canonical(const Eigen::Matrix3d& image_conic);
}
