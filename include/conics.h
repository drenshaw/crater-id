#pragma once

#include <iostream>
#include <tuple>
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include "/usr/local/include/vtk-9.2/vtk3DS.h"

#include "vector_math.h"

class Conic {
    public:
        Conic(const double=0, const double=0, const double=0, const double=0, const double=0);
        Conic(const std::tuple<double, double, double, double, double>&);
        Conic(const std::vector<double>&);
        void setGeometricParameters(const std::vector<double>&);
        void setImplicitParameters(const std::vector<double>& impl_params);
        void setLocus(const Eigen::Matrix3d& locus);
        void NormalizeImplicitParameters(std::vector<double>&);
        Eigen::Vector2d getCenter();
        Eigen::Vector2d getSemiAxes();
        std::vector<double> getGeom();
        std::vector<double> getImplicit();
        Eigen::Matrix3d getLocus();
        Eigen::Matrix3d getEnvelope();
        // Eigen::Matrix3d getMatrixAdjugate(const Eigen::Matrix3d&);
        std::vector<double> Locus2Implicit(const Eigen::Matrix3d&);
        std::vector<double> Implicit2Geom(const std::vector<double>&);
        Eigen::Matrix3d Geom2Locus();
        Eigen::Matrix3d Implicit2Locus(const std::vector<double>&);
        std::vector<double> Locus2Geom(const Eigen::Matrix3d&);
        std::vector<double> Geom2Implicit();
        bool ConicIntersectionLines(const Eigen::Matrix3d&, 
                                    std::tuple<Eigen::Vector3d, Eigen::Vector3d>&);
        bool ConicIntersectionLines(Conic&,
                                    std::tuple<Eigen::Vector3d, Eigen::Vector3d>&);
     

    private:
        double semimajor_axis_;
        double semiminor_axis_;
        double x_center_;
        double y_center_;
        double angle_;
        Eigen::Matrix3d locus_;
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

Eigen::Matrix3d getMatrixAdjugate(const Eigen::Matrix3d&);
bool IntersectionLines(const Eigen::Matrix3d&, 
                       const Eigen::Matrix3d&,
                       std::tuple<Eigen::Vector3d, Eigen::Vector3d>&);
bool ChooseIntersection(const std::tuple<Eigen::Vector3d, Eigen::Vector3d>&, 
                        const Eigen::Vector2d&, const Eigen::Vector2d&, Eigen::Vector3d&);
bool IntersectConics(const Eigen::Matrix3d&, 
                     const Eigen::Matrix3d&, 
                     const double,
                     std::tuple<Eigen::Vector3d, Eigen::Vector3d>&);
bool computeInvariant(const Eigen::Vector3d&, 
                      const Eigen::Vector3d&, 
                      const Eigen::Matrix3d&,
                      double&);
bool computeCraterTriadInvariants(Conic& A, Conic& B, Conic& C,
                                  std::vector<double>& invariants);
void cylinder();