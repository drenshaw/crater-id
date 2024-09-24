#ifndef QUADRICS_H
#define QUADRICS_H

#include <iostream>
#include <eigen3/Eigen/Dense>

class Quadric {
  public:
    Quadric(const Eigen::Vector3d& position, const double radius, const Eigen::Vector3d& surface_normal, const std::string id="");
    Quadric(const Eigen::Vector3d& position, const double radius, const std::string id="");
    Quadric(const double lat, const double lon, const double radius, const std::string id="");
    Quadric(const std::string id, const Eigen::Vector3d&, const Eigen::Vector3d&, const double);
    Quadric(const std::string id, const Eigen::Vector3d&, const double);
    Quadric(const std::string id, const double, const double, const double);
    Quadric(const std::string id, const Eigen::Vector3d& position, const double radius, const Eigen::Vector3d& surface_normal);
    Eigen::Matrix3d getQuadricTransformationMatrix() const;
    Eigen::Matrix4d getLocus() const;
    Eigen::Vector3d getLocation() const;
    void getLocation(Eigen::Vector3d& location) const;
    Eigen::Vector3d getNormal() const;
    void getNormal(Eigen::Vector3d& surface_normal) const;
    Eigen::Vector4d getPlane() const;
    void getPlane(Eigen::Vector4d& surface_normal) const;
    double getAngleBetweenQuadrics(const Quadric& other_quadric) const;
    Eigen::Vector3d getAxisNormalToQuadrics(const Quadric& other_quadric) const;

  private:
    Eigen::Matrix4d generateQuadricLocus() const ;    
    void loadQuadric(const Eigen::Matrix3d& conic_locus);
    // All positions or directions are given in the Moon-centered frame
    // aka, the 'selenographic' frame
    std::string id_;
    double radius_;
    Eigen::Vector3d surface_point_;
    Eigen::Vector3d surface_normal_;
    Eigen::Vector4d plane_;
    
    friend std::ostream& operator<<(std::ostream& os, const Quadric&);
};


double calculateCraterRimFromRadius(const double radius);
Eigen::Matrix4d GenerateQuadricLocusFromRadiusNormal(const Eigen::Vector3d& position, const double radius);
Eigen::Matrix4d ConicEnvelopeToQuadricEnvelope(const Eigen::Matrix3d& conic_envelope, 
                                               const Eigen::MatrixXd& h_k);
Eigen::Vector4d SurfacePointToPlane(const Eigen::Matrix3d& T_e2m, 
                                    const Eigen::Vector3d& surface_point);
void GenerateQuadricFromRadiusNormal();
Eigen::MatrixXd transformSelenographicToCraterFrame(const Eigen::Vector3d&, 
                                                    const Eigen::Matrix3d& T_e2m);

#endif