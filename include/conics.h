#ifndef CONICS_H
#define CONICS_H

#include <iostream>
#include <stdexcept>
#include <tuple>
#include <array>
#include <math.h>
#include <Eigen/Dense>
#include <opencv2/core/core.hpp>

#include "vector_math.h"

#define GEOMETRIC_PARAM 5
#define IMPLICIT_PARAM 6
#define CONIC_DIM 3

template <typename T>
bool almost_equal(const T a, const T b) {
  return std::abs(a-b)<EPS;
}

class Conic {
  public:
    Conic(const double=0, const double=0, const double=0, 
          const double=0, const double=0);
    Conic(const std::array<double, GEOMETRIC_PARAM>&);
    bool operator==(Conic& other_conic);
    bool operator!=(const Conic& other_conic);
    void SetGeometricParameters(const std::array<double, GEOMETRIC_PARAM>&);
    void SetGeometricParameters(const std::vector<double>&);
    void SetGeometricParameters(const double, 
                                const double, 
                                const double, 
                                const double, 
                                const double);
    void SetImplicitParameters(const std::array<double, IMPLICIT_PARAM>&);
    void SetLocus(const Eigen::Matrix3d& locus);
    void NormalizeImplicitParameters(std::vector<double>&);
    void NormalizeImplicitParameters(std::array<double, IMPLICIT_PARAM>&);
    Eigen::Vector2d GetCenter();
    double GetCenterX();
    double GetCenterY();
    void GetCenter(Eigen::Vector2d& center);
    void GetCenter(cv::Point& center);
    double GetSemiMajorAxis();
    double GetSemiMinorAxis();
    Eigen::Vector2d GetSemiAxes();
    void GetSemiAxes(Eigen::Vector2d& semiaxes);
    void GetSemiAxes(cv::Point& semiaxes);
    double GetAngle();
    int GetID();
    std::array<double, GEOMETRIC_PARAM> GetGeom();
    std::array<double, IMPLICIT_PARAM> GetImplicit();
    Eigen::Matrix3d GetLocus();
    Eigen::Matrix3d GetEnvelope();
    // Eigen::Matrix3d getMatrixAdjugate(const Eigen::Matrix3d&);
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
    double semimajor_axis_;
    double semiminor_axis_;
    double x_center_;
    double y_center_;
    double angle_;
    unsigned int id_;
    void setID();
    void MakeConic( const double semimajor_axis, 
                const double semiminor_axis, 
                const double x_center, 
                const double y_center, 
                const double angle);
    void MakeConic(const std::array<double, GEOMETRIC_PARAM>& geom_arr);
    void MakeConic(const std::vector<double>& geom_vec);

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

Eigen::Vector3d GetNorthPoleUnitVector();
void GetNorthPoleUnitVector(Eigen::Vector3d&);
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
bool computeCraterTriadInvariants(const Conic&, const Conic&, const Conic&,
                                  std::array<double, CONIC_DIM>&);
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
#endif