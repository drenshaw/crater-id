#ifndef QUADRICS_H
#define QUADRICS_H

#include <iostream>
#include "conics.h"
#include "vector_math.h"

class Quadric {
  public:
    Quadric(const Eigen::Vector3d& position, const double radius, const Eigen::Vector3d& surface_normal, const std::string id="");
    Quadric(const Eigen::Vector3d& position, const double radius, const std::string id="");
    Quadric(const double lat, const double lon, const double radius, const std::string id="");
    Quadric(const std::string id, const Eigen::Vector3d&, const Eigen::Vector3d&, const double);
    Quadric(const std::string id, const Eigen::Vector3d&, const double);
    Quadric(const std::string id, const double, const double, const double);
    Quadric(const std::string id, const Eigen::Vector3d& position, const double radius, const Eigen::Vector3d& surface_normal);
    Eigen::Matrix3d GetQuadricTransformationMatrix();
    Eigen::Matrix4d GetLocus();

  private:
    Eigen::Matrix4d GenerateQuadricLocus();
    void MakeQuadric(const Eigen::Vector3d& position, const double radius, const Eigen::Vector3d& surface_normal, const std::string id);
    void MakeQuadric(const Eigen::Vector3d& position, const double radius, const std::string id);
    void MakeQuadric(const double lat, const double lon, const double radius, const std::string id);
    void LoadQuadric(const Eigen::Matrix3d& conic_locus);
    // All positions or directions are given in the Moon-centered frame
    // aka, the 'selenographic' frame
    std::string id_;
    double radius_;
    Eigen::Vector3d surface_point_;
    Eigen::Vector3d surface_normal_;
    Eigen::Vector3d plane_normal_;
    Eigen::Vector4d plane_;
    Eigen::Matrix3d T_e2m_;
    Eigen::Matrix4d locus_;
    Eigen::Matrix4d envelope_;
    
    friend std::ostream& operator<<(std::ostream& os, const Quadric&);
};

Eigen::Matrix4d GenerateQuadricFromRadiusNormal(const Eigen::Vector3d& position, const double radius);
Eigen::Matrix4d ConicEnvelopeToQuadricEnvelope(const Eigen::Matrix3d& conic_envelope, 
                                               const Eigen::MatrixXd& h_k);
std::tuple<Eigen::Vector4d, Eigen::Vector3d> 
SurfacePointToPlane(const Eigen::Matrix3d& T_e2m, 
                    const Eigen::Vector3d& surface_point);

#endif