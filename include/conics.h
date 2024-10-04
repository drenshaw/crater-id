#pragma once

// #include <iostream>
// // #include <stdexcept>
// #include <tuple>
// #include <array>
// #include <math.h>
// #include <eigen3/Eigen/Dense>
#include <opencv2/core/core.hpp>

#include "math_utils.h"

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
    Conic(const double semimajor_axis=1, const double semiminor_axis=1, 
          const double x_center=0, const double y_center=0, const double angle=0);
    Conic(const std::array<double, GEOMETRIC_PARAM>&);
    Conic(const std::vector<double>&);
    Conic(const Eigen::Matrix3d& locus);
    bool operator==(const Conic& other_conic) const;
    bool operator!=(const Conic& other_conic) const;
    bool operator==(const Conic* other_conic) const;
    bool operator!=(const Conic* other_conic) const;
    void setGeometricParameters(const std::array<double, GEOMETRIC_PARAM>&);
    void setGeometricParameters(const std::vector<double>&);
    void setGeometricParameters(const double, 
                                const double, 
                                const double, 
                                const double, 
                                const double);
    void setSemimajorAxis(const double semimajor_axis);
    void setSemiminorAxis(const double semiminor_axis);
    void setCenterX(const double x_center);
    void setCenterY(const double y_center);
    void setAngle(const double angle);                                
    void setImplicitParameters(const std::array<double, IMPLICIT_PARAM>&);
    void setLocus(const Eigen::Matrix3d& locus);
    Eigen::Vector2d getCenter() const;
    double getCenterX() const;
    double getCenterY() const;
    void getCenter(Eigen::Vector2d& center) const;
    void getCenter(cv::Point& center) const;
    double getSemiMajorAxis() const;
    double getSemiMinorAxis() const;
    cv::Size getSemiAxes() const;
    void getSemiAxes(Eigen::Vector2d& semiaxes) const;
    void getSemiAxes(cv::Point& semiaxes) const;
    void getSemiAxes(cv::Size& semiaxes) const;
    cv::Size getSize() const;
    void getSize(cv::Size& semiaxes) const;
    double getAngle() const;
    int getID() const;
    void normalizeImplicitParams();
    std::array<double, GEOMETRIC_PARAM> getGeom() const ;
    std::array<double, IMPLICIT_PARAM> getImplicit() const ;
    Eigen::Matrix3d getLocus() const ;
    Eigen::Matrix3d getEnvelope() const ;
    Eigen::Matrix3d toLocus() const ;
    std::array<double, GEOMETRIC_PARAM> fromLocus(const Eigen::Matrix3d&) const ;
    std::array<double, IMPLICIT_PARAM> toImplicit() const ;
    bool intersectsConic(const Eigen::Matrix3d&, 
                                std::tuple<Eigen::Vector3d, Eigen::Vector3d>&) const ;
    bool intersectsConicLines(const Conic&,
                                std::tuple<Eigen::Vector3d, Eigen::Vector3d>&) const ;
    bool chooseIntersection(const Conic& other, Eigen::Vector3d&) const ;
    
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


/*********************************************************/
/***********************Conic Utils***********************/
/*********************************************************/

// Convert to/from representations
std::array<double, IMPLICIT_PARAM> locus2Implicit(const Eigen::Matrix3d& locus);
std::array<double, GEOMETRIC_PARAM> implicit2Geom(const std::array<double, IMPLICIT_PARAM>& impl_params);
std::array<double, IMPLICIT_PARAM> geom2Implicit( const double semimajor_axis, 
                                                  const double semiminor_axis, 
                                                  const double x_center, 
                                                  const double y_center, 
                                                  const double angle);

std::array<double, GEOMETRIC_PARAM> locus2Geom(const Eigen::Matrix3d& locus);
Eigen::Matrix3d implicit2Locus(const std::array<double, IMPLICIT_PARAM>& impl_params);

 // General Conic Utils
void normalizeImplicitParameters(std::array<double, IMPLICIT_PARAM>& impl_params) ;
void normalizeImplicitParameters(std::vector<double>& impl_params) ;


/*********************************************************/
/***********************INVARIANTS************************/
/*********************************************************/
namespace invariants {
// TODO: This function probably doesn't need to be public
Eigen::Matrix3d stackVectorsInto3x3Matrix(const Eigen::Vector3d& A, 
                                          const Eigen::Vector3d& B, 
                                          const Eigen::Vector3d& C);
double crossRatio(const Eigen::Vector3d& ref, 
                  const Eigen::Vector3d& A, 
                  const Eigen::Vector3d& B, 
                  const Eigen::Vector3d& C, 
                  const Eigen::Vector3d& D);
bool intersectionLines(const Eigen::Matrix3d&, 
                       const Eigen::Matrix3d&,
                       std::tuple<Eigen::Vector3d, Eigen::Vector3d>&);
bool chooseIntersection(const std::tuple<Eigen::Vector3d, 
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
                                  std::array<double, NONCOPLANAR_INVARIANTS>&);
} // namespace
