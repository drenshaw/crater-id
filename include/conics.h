#pragma once

#include <iostream>
#include <tuple>
#include <vector>
#include <math.h>
// #include <cmath>
#include <Eigen/Dense>
// #include <Eigen/Eigenvalues>

#include "vector_math.h"

// template <typename T>
class Conic {
    public:
        Conic(const float=0, const float=0, const float=0, const float=0, const float=0);
        Conic(const std::tuple<float, float, float, float, float>&);
        Conic(const std::vector<float>&);
        void setGeometricParameters(const std::vector<float>&);
        void setImplicitParameters(const std::vector<float>& impl_params);
        void setLocus(const Eigen::Matrix3f& locus);
        void NormalizeImplicitParameters(std::vector<float>&);
        Eigen::Vector2f getCenter();
        Eigen::Vector2f getSemiAxes();
        std::vector<float> getGeom();
        std::vector<float> getImplicit();
        Eigen::Matrix3f getLocus();
        Eigen::Matrix3f getEnvelope();
        // Eigen::Matrix3f getMatrixAdjugate(const Eigen::Matrix3f&);
        std::vector<float> Locus2Implicit(const Eigen::Matrix3f&);
        std::vector<float> Implicit2Geom(const std::vector<float>&);
        Eigen::Matrix3f Geom2Locus();
        Eigen::Matrix3f Implicit2Locus(const std::vector<float>&);
        std::vector<float> Locus2Geom(const Eigen::Matrix3f&);
        std::vector<float> Geom2Implicit();
        bool ConicIntersectionLines(const Eigen::Matrix3f&, 
                                    std::tuple<Eigen::Vector3f, Eigen::Vector3f>&);
        bool ConicIntersectionLines(Conic&,
                                    std::tuple<Eigen::Vector3f, Eigen::Vector3f>&);
     

    private:
        float semimajor_axis_;
        float semiminor_axis_;
        float x_center_;
        float y_center_;
        float angle_;
        // std::tuple<float, float, float, float, float> geom;
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

Eigen::Matrix3f getMatrixAdjugate(const Eigen::Matrix3f&);
bool IntersectionLines(const Eigen::Matrix3f&, 
                       const Eigen::Matrix3f&,
                       std::tuple<Eigen::Vector3f, Eigen::Vector3f>&);
bool ChooseIntersection(const std::tuple<Eigen::Vector3f, Eigen::Vector3f>&, 
                        const Eigen::Vector2f&, const Eigen::Vector2f&, Eigen::Vector3f&);
bool IntersectConics(const Eigen::Matrix3f&, 
                     const Eigen::Matrix3f&, 
                     const float,
                     std::tuple<Eigen::Vector3f, Eigen::Vector3f>&);
bool computeInvariant(const Eigen::Vector3f&, 
                      const Eigen::Vector3f&, 
                      const Eigen::Matrix3f&,
                      float&);
