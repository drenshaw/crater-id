#ifndef QUADRICS_H
#define QUADRICS_H

#include <iostream>
#include <eigen3/Eigen/Geometry>

class Quadric {
  public:
    Quadric(const Eigen::Vector3d& position, const double radius, const Eigen::Vector3d& surface_normal, const std::string id="");
    Quadric(const Eigen::Vector3d& position, const double radius, const std::string id="");
    Quadric(const double lat, const double lon, const double radius, const std::string id="");
    Quadric(const std::string id, const Eigen::Vector3d&, const Eigen::Vector3d&, const double);
    Quadric(const std::string id, const Eigen::Vector3d&, const double);
    Quadric(const std::string id, const double, const double, const double);
    Quadric(const std::string id, const Eigen::Vector3d& position, const double radius, const Eigen::Vector3d& surface_normal);
    Quadric(const Eigen::Vector3d& pt1, const Eigen::Vector3d& pt2, const Eigen::Vector3d& pt3, std::string id="");
    
    bool operator==(const Quadric& other_quadric) const;
    bool operator!=(const Quadric& other_quadric) const;
    bool operator==(const Quadric* other_quadric) const;
    bool operator!=(const Quadric* other_quadric) const;
    
    Eigen::Matrix3d getQuadricTransformationMatrix() const;
    Eigen::Matrix4d getLocus() const;
    Eigen::Vector3d getLocation() const;
    void getLocation(Eigen::Vector3d& location) const;
    Eigen::Vector3d getNormal() const;
    void getNormal(Eigen::Vector3d& surface_normal) const;
    Eigen::Hyperplane<double, 3> getPlane() const;
    void getPlane(Eigen::Hyperplane<double, 3>& hyperplane) const;
    double getAngleBetweenQuadrics(const Quadric& other_quadric) const;
    Eigen::Vector3d getAxisNormalToQuadrics(const Quadric& other_quadric) const;
    double getRadius() const;
    std::string getID() const;

  private:
    Eigen::Matrix4d generateQuadricLocus() const ;    
    void loadQuadric(const Eigen::Matrix3d& conic_locus);
    // All positions or directions are given in the Moon-centered frame
    // aka, the 'selenographic' frame
    std::string id_;
    double radius_;
    Eigen::Vector3d surface_point_;
    Eigen::Hyperplane<double, 3> plane_;
    
    friend std::ostream& operator<<(std::ostream& os, const Quadric&);
};


double calculateCraterRimFromRadius(const double radius);
Eigen::Matrix4d GenerateQuadricLocusFromRadiusNormal(const Eigen::Vector3d& position, const double radius);
Eigen::Matrix4d ConicEnvelopeToQuadricEnvelope(const Eigen::Matrix3d& conic_envelope, 
                                               const Eigen::MatrixXd& h_k);
Eigen::Hyperplane<double, 3> SurfacePointToPlane(const Eigen::Matrix3d& T_e2m, 
                                                      const Eigen::Vector3d& surface_point);
void GenerateQuadricFromRadiusNormal();
Eigen::MatrixXd transformSelenographicToCraterFrame(const Eigen::Vector3d&, 
                                                    const Eigen::Matrix3d& T_e2m);

#endif