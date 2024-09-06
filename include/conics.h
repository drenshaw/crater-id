#ifndef CONICS_H
#define CONICS_H

#include <iostream>
#include <stdexcept>
#include <tuple>
#include <array>
#include <math.h>
#include <eigen3/Eigen/Dense>
#include <opencv2/core/core.hpp>

#include "vector_math.h"

#define GEOMETRIC_PARAM 5
#define IMPLICIT_PARAM 6
#define CONIC_DIM 3
#define NONCOPLANAR_INVARIANTS 3
#define    COPLANAR_INVARIANTS 7

template <typename T>
bool almost_equal(const T a, const T b) {
  return std::abs(a-b)<EPS;
}

class Conic {
  public:
  // TODO: update constructor with modern C++ member initializers
    Conic(const double=1, const double=1, const double=0, 
          const double=0, const double=0);
    Conic(const std::array<double, GEOMETRIC_PARAM>&);
    Conic(const std::vector<double>&);
    bool operator==(const Conic& other_conic) const;
    bool operator!=(const Conic& other_conic) const;
    void SetGeometricParameters(const std::array<double, GEOMETRIC_PARAM>&);
    void SetGeometricParameters(const std::vector<double>&);
    void SetGeometricParameters(const double, 
                                const double, 
                                const double, 
                                const double, 
                                const double);
    void SetSemimajorAxis(const double semimajor_axis);
    void SetSemiminorAxis(const double semiminor_axis);
    void SetCenterX(const double x_center);
    void SetCenterY(const double y_center);
    void SetAngle(const double angle);                                
    void SetImplicitParameters(const std::array<double, IMPLICIT_PARAM>&);
    void SetLocus(const Eigen::Matrix3d& locus);
    void NormalizeImplicitParameters(std::array<double, IMPLICIT_PARAM>&);
    void NormalizeImplicitParameters(std::vector<double>&);
    Eigen::Vector2d GetCenter() const;
    double GetCenterX() const;
    double GetCenterY() const;
    void GetCenter(Eigen::Vector2d& center) const;
    void GetCenter(cv::Point& center) const;
    double GetSemiMajorAxis() const;
    double GetSemiMinorAxis() const;
    Eigen::Vector2d GetSemiAxes() const;
    void GetSemiAxes(Eigen::Vector2d& semiaxes) const;
    void GetSemiAxes(cv::Point& semiaxes) const;
    double GetAngle() const;
    int GetID() const;
    std::array<double, GEOMETRIC_PARAM> GetGeom();
    std::array<double, IMPLICIT_PARAM> GetImplicit();
    Eigen::Matrix3d GetLocus();
    Eigen::Matrix3d GetEnvelope();
    // Eigen::Matrix3d getAdjugateMatrix(const Eigen::Matrix3d&);
    std::array<double, IMPLICIT_PARAM> Locus2Implicit(const Eigen::Matrix3d&);
    std::array<double, GEOMETRIC_PARAM> Implicit2Geom(const std::array<double, IMPLICIT_PARAM>&);
    Eigen::Matrix3d Geom2Locus();
    Eigen::Matrix3d Implicit2Locus(const std::array<double, IMPLICIT_PARAM>&);
    Eigen::Matrix3d Implicit2Locus();
    std::array<double, GEOMETRIC_PARAM> Locus2Geom(const Eigen::Matrix3d&);
    std::array<double, IMPLICIT_PARAM> Geom2Implicit();
    bool ConicIntersectionLines(const Eigen::Matrix3d&, 
                                std::tuple<Eigen::Vector3d, Eigen::Vector3d>&);
    bool ConicIntersectionLines(Conic&,
                                std::tuple<Eigen::Vector3d, Eigen::Vector3d>&);
    
  protected:
    static int next_id;

  private:
    double semimajor_axis_ = 1;
    double semiminor_axis_ = 1;
    double x_center_ = 0;
    double y_center_ = 0;
    double angle_ = 0;
    unsigned int id_;
    void setID();

    friend std::ostream& operator<<(std::ostream& os, const Conic&);
};

class ConicImplicit {
  ConicImplicit();
};

class ConicGeometry {
  ConicGeometry();
};

class ConicMatrix {
  ConicMatrix();
}; 

bool IntersectionLines(const Eigen::Matrix3d&, 
                       const Eigen::Matrix3d&,
                       std::tuple<Eigen::Vector3d, Eigen::Vector3d>&);
bool ChooseIntersection(const std::tuple<Eigen::Vector3d, 
                        Eigen::Vector3d>&, 
                        const Eigen::Vector2d&, 
                        const Eigen::Vector2d&, 
                        Eigen::Vector3d&);
bool IntersectConics(const Eigen::Matrix3d&, 
                     const Eigen::Matrix3d&, 
                     const double,
                     std::tuple<Eigen::Vector3d, Eigen::Vector3d>&);
bool computeInvariant(const Eigen::Vector3d&, 
                      const Eigen::Vector3d&, 
                      const Eigen::Matrix3d&,
                      double&);
bool computeCraterTriadInvariants(Conic&, Conic&, Conic&,
                                  std::array<double, NONCOPLANAR_INVARIANTS>&);
Eigen::Matrix3d getENUFrame(const Eigen::Vector3d&);
void GenerateQuadricFromRadiusNormal();
Eigen::MatrixXd transformSelenographicToCraterFrame(const Eigen::Vector3d&, 
                                                    const Eigen::Matrix3d& T_e2m);
Eigen::Matrix3d pointCameraInDirection(const Eigen::Vector3d& camera_position, 
                                       const Eigen::Vector3d& desired_location);
Eigen::Quaterniond eulerToQuaternion(const double roll, const double pitch, const double yaw);
void eulerToQuaternion(const double roll,
                       const double pitch,
                       const double yaw,
                       Eigen::Matrix3d& dcm);    

void convertEigenVectorToVector(const Eigen::Vector3d& eig, std::array<double, CONIC_DIM>& arr);
void convertEigenVectorToVector(const Eigen::Vector3d& eig, std::vector<double>& vec);
bool vectorContainsNaN(const Eigen::Vector3d& eV);
template <typename T, size_t SIZE>
bool vectorContainsNaN(const std::array<T, SIZE>& vec) {
  return std::any_of(vec.begin(), vec.end(), [](T i){return std::isnan(i);});
}

#endif