#include "quadrics.h"
#include "io.h"
#include "vector_math.h"
#include "conics.h"

// This is the base constructor
Quadric::Quadric( const Eigen::Vector3d& position, 
                  const double radius, 
                  const Eigen::Vector3d& surface_normal, 
                  const std::string id) : 
                    id_(id),
                    radius_(radius),
                    surface_point_{position},
                    surface_normal_(surface_normal.normalized()) {
  // TODO: Uses only Moon radius; change if using another body (e.g., Mars)
  if(radius > R_MOON) {
    throw std::runtime_error("The crater radius is larger than the Moon.");
  }
  if(std::abs(surface_normal.norm()) < 1e-8) {
    throw std::runtime_error("The surface normal is not defined or ==0.");
  }
  if(radius == 0) {
    throw std::runtime_error("The crater radius is zero.");
  }
  Eigen::Matrix3d T_e2m_ = GetQuadricTransformationMatrix();
  plane_ = SurfacePointToPlane(T_e2m_, surface_point_);
}

// These are delegating constructors
Quadric::Quadric(const Eigen::Vector3d& position, const double radius, const std::string id) : 
  Quadric(position, radius, position.normalized(), id) {}
Quadric::Quadric(const std::string id, const Eigen::Vector3d& position, const double radius) :
  Quadric(position, radius, id) {}
Quadric::Quadric(const double lat, const double lon, const double radius, const std::string id) :
  Quadric(calculateCraterRimFromRadius(radius)*latlon2bearing(lat, lon), radius, id) {}
Quadric::Quadric(const std::string id, const double lat, const double lon, const double radius) :
  Quadric(calculateCraterRimFromRadius(radius)*latlon2bearing(lat, lon), radius, id) {}
Quadric::Quadric(const std::string id, const Eigen::Vector3d& position, const double radius, const Eigen::Vector3d& surface_normal) :
  Quadric(position, radius, surface_normal, id) {}
        
Eigen::Matrix4d Quadric::GenerateQuadricLocus() {
  return GenerateQuadricFromRadiusNormal(surface_point_, radius_);
}

Eigen::Matrix3d Quadric::GetQuadricTransformationMatrix() {
  return getENUFrame(surface_normal_);
}

Eigen::Matrix4d Quadric::GetLocus() {
  return GenerateQuadricLocus();
}

std::ostream& operator<<(std::ostream& os, const Quadric& quad) {
  return os 
    << "Quadric -> \tID: '" << quad.id_ << "'"
    << "\n\t\tLocation: (" << quad.surface_point_.transpose() << ") | "
    << "\tRadius: " << quad.radius_
    << "\n\t\tPlane: [" << quad.plane_.transpose() << "] ";
}

double calculateCraterRimFromRadius(const double radius) {
  return sqrt(pow(R_MOON, 2) - pow(radius, 2));
}

Eigen::Vector4d SurfacePointToPlane(const Eigen::Matrix3d& T_e2m, 
                                    const Eigen::Vector3d& surface_point) {
  Eigen::Vector3d u_north_pole = getNorthPoleUnitVector();
  Eigen::Vector3d plane_normal = T_e2m * u_north_pole;
  double rho = surface_point.dot(plane_normal);
  Eigen::Vector4d plane;
  plane << plane_normal, -rho;
  return plane;
}

Eigen::Matrix4d GenerateQuadricFromRadiusNormal(const Eigen::Vector3d& position, const double radius) {
  Conic conic(radius, radius, 0, 0, 0);
  Eigen::Matrix3d conic_envelope = conic.GetEnvelope();

  Eigen::Matrix3d T_enu_to_ref = getENUFrame(position);
  Eigen::MatrixXd h_k = transformSelenographicToCraterFrame(position, T_enu_to_ref);
  // eq 40 of Christian, Derksen, and Watkins [2020]
  Eigen::Matrix4d quadric_envelope = ConicEnvelopeToQuadricEnvelope(conic_envelope, h_k);
  Eigen::Matrix4d quadric_locus = getAdjugateMatrix(quadric_envelope);
  // bool success = normalizeDeterminant(quadric_locus);
  return quadric_locus;
}

Eigen::Matrix4d ConicEnvelopeToQuadricEnvelope(const Eigen::Matrix3d& conic_envelope, const Eigen::MatrixXd& h_k) {
  // transformation from crater plane to selenographic plane
  // convert to homogeneous coordinates; 4x3
  return h_k * conic_envelope * h_k.transpose();
}
