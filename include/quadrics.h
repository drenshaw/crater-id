#pragma once

#include "conics.h"
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
    Eigen::Matrix4d getEnvelope() const;
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
    Eigen::Matrix3d projectToImageEnvelope(const Eigen::MatrixXd& proj_mtx) const;
    Eigen::Matrix3d projectToImageLocus(const Eigen::MatrixXd& proj_mtx) const;
    Conic projectToImage(const Eigen::MatrixXd& proj_mtx) const;
    Eigen::Matrix3d projectToPlaneEnvelope(const Eigen::MatrixXd& extrinsic_mtx) const;
    Eigen::Matrix3d projectToPlaneLocus(const Eigen::MatrixXd& extrinsic_mtx) const;
    Conic projectToImagePlane(const Eigen::MatrixXd& extrinsic_mtx) const;

  private:
    Eigen::Matrix4d generateQuadricLocusFromPointRadius() const;
    Eigen::Matrix4d generateQuadricEnvelopeFromPointRadius() const;
    // All positions or directions are given in the Moon-centered frame
    // aka, the 'selenographic' frame
    std::string id_;
    double radius_;
    Eigen::Vector3d center_;
    Eigen::Hyperplane<double, 3> plane_;
    
    friend std::ostream& operator<<(std::ostream& os, const Quadric&);
};

bool isSamePlane(const Eigen::Hyperplane<double, 3>& p1, const Eigen::Hyperplane<double, 3>& p2, const double thresh=1e-3);
bool isSamePlane(const Quadric& quad1, const Quadric& quad2, const double thresh=1e-3);
Eigen::Matrix4d GenerateQuadricLocus(const Eigen::Vector3d& position, const double radius);
Eigen::Matrix4d GenerateQuadricEnvelope(const Eigen::Vector3d& position, const double radius);
Eigen::Matrix4d ConicEnvelopeToQuadricEnvelope(const Eigen::Matrix3d& conic_envelope, 
                                               const Eigen::MatrixXd& h_k);
Eigen::Hyperplane<double, 3> SurfacePointToPlane(const Eigen::Matrix3d& T_e2m, 
                                                      const Eigen::Vector3d& surface_point);
Eigen::MatrixXd transformSelenographicToCraterFrame(const Eigen::Vector3d&, 
                                                    const Eigen::Matrix3d& T_e2m);
Eigen::Matrix3d getAttitudeTransformBetweenPoints(const Eigen::Vector3d& camera_position, 
                                       const Eigen::Vector3d& desired_location);


Eigen::Matrix4d makeSphere(const double radius);
Eigen::Matrix4d makeEllipsoid(const Eigen::Vector3d& radii);                                       
