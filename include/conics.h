#ifndef CONICS_H
#define CONICS_H

#include <iostream>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <opencv2/core/core.hpp>

#include "vector_math.h"

#define GEOMETRIC_PARAM 5
#define IMPLICIT_PARAM 6

template <typename T>
bool almost_equal(const T a, const T b) {
  return std::abs(a-b)<EPS;
}

class Conic {
  public:
    typedef std::tuple<double, double, double, double, double, double> tuple_impl;
    // typedef std::tuple<double, double, double, double, double, double> tuple_impl;
    Conic(const double=0, const double=0, const double=0, 
          const double=0, const double=0);
    Conic(const std::tuple<double, double, double, double, double>&);
    Conic(const std::vector<double>&);
    bool operator==(const Conic&);
    bool operator!=(const Conic& other_conic);
    void setGeometricParameters(const std::vector<double>&);
    void setImplicitParameters(const std::vector<double>& impl_params);
    void setLocus(const Eigen::Matrix3d& locus);
    void NormalizeImplicitParameters(std::vector<double>&);
    Eigen::Vector2d getCenter();
    void getCenter(Eigen::Vector2d& center);
    void getCenter(cv::Point& center);
    Eigen::Vector2d getSemiAxes();
    void getSemiAxes(Eigen::Vector2d& semiaxes);
    void getSemiAxes(cv::Point& semiaxes);
    double getAngle();
    int getID();
    std::vector<double> getGeom();
    std::vector<double> getImplicit();
    Eigen::Matrix3d getLocus();
    Eigen::Matrix3d getEnvelope();
    // Eigen::Matrix3d getMatrixAdjugate(const Eigen::Matrix3d&);
    std::vector<double> Locus2Implicit(const Eigen::Matrix3d&);
    std::vector<double> Implicit2Geom(const std::vector<double>&);
    Eigen::Matrix3d Geom2Locus();
    Eigen::Matrix3d Implicit2Locus(const std::vector<double>&);
    Eigen::Matrix3d Implicit2Locus();
    std::vector<double> Locus2Geom(const Eigen::Matrix3d&);
    std::vector<double> Geom2Implicit();
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
    void conic(const double, const double, const double, const double, const double);
    // Eigen::Matrix3d locus_;
    // std::vector<double> implicit_;
    // std::tuple<double, double, double, double, double> geom;
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

Eigen::Vector3d getNorthPoleUnitVector();
void getNorthPoleUnitVector(Eigen::Vector3d&);
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
                                  std::vector<double>&);
Eigen::Matrix3d getENUFrame(const Eigen::Vector3d&);
void generateQuadricFromRadiusNormal();
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