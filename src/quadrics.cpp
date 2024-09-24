#include "quadrics.h"
#include "io.h"
#include "math_utils.h"
#include "conics.h"

#include <iomanip>

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
  if(surface_normal.norm() < 1e-8) {
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
        
Eigen::Matrix4d Quadric::GenerateQuadricLocus() const {
  return GenerateQuadricLocusFromRadiusNormal(surface_point_, radius_);
}

Eigen::Matrix3d Quadric::GetQuadricTransformationMatrix() const {
  return getENUFrame(surface_normal_);
}

Eigen::Matrix4d Quadric::GetLocus() const {
  return GenerateQuadricLocus();
}

Eigen::Vector3d Quadric::GetLocation() const {
  return surface_point_;
}

void Quadric::GetLocation(Eigen::Vector3d& location) const {
  location = surface_point_;
}

Eigen::Vector3d Quadric::GetNormal() const {
  return surface_normal_;
}

void Quadric::GetNormal(Eigen::Vector3d& surface_normal) const {
  surface_normal = surface_normal_;
}

Eigen::Vector4d Quadric::GetPlane() const {
  return plane_;
}

void Quadric::GetPlane(Eigen::Vector4d& plane) const {
  plane = plane_;
}


std::ostream& operator<<(std::ostream& os, const Quadric& quad) {
  std::streamsize ss = std::cout.precision();
  std::streamsize sw = std::cout.width();
  return os 
    << "Quadric -> \tID: '" << quad.id_ << "'"
    << std::fixed << std::setw(8) << std::setprecision(1)
    << "\n\t\tLocation: (" << quad.surface_point_.transpose() << ") km | "
    << "\tRadius: " << quad.radius_ << "km"
    << "\n\t\tPlane: [" << quad.plane_.transpose() << " ] "
    << std::setprecision(ss) << std::setw(sw) << std::defaultfloat
    ;
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

Eigen::Matrix4d GenerateQuadricLocusFromRadiusNormal(const Eigen::Vector3d& position, const double radius) {
  Conic conic(radius, radius, 0, 0, 0);
  Eigen::Matrix3d conic_envelope = conic.GetEnvelope();

  Eigen::Matrix3d T_enu_to_ref = getENUFrame(position);
  std::cout << T_enu_to_ref << std::endl;
  Eigen::MatrixXd h_k = transformSelenographicToCraterFrame(position, T_enu_to_ref);
  // eq 40 of Christian, Derksen, and Watkins [2020]
  Eigen::Matrix4d quadric_envelope = ConicEnvelopeToQuadricEnvelope(conic_envelope, h_k);
  Eigen::Matrix4d quadric_locus = getAdjugateMatrix(quadric_envelope);
  // bool success = normalizeDeterminant(quadric_locus);
  return quadric_locus/quadric_locus(3,3);
}

Eigen::Matrix4d ConicEnvelopeToQuadricEnvelope(const Eigen::Matrix3d& conic_envelope, const Eigen::MatrixXd& h_k) {
  // transformation from crater plane to selenographic plane
  // convert to homogeneous coordinates; 4x3
  return h_k * conic_envelope * h_k.transpose();
}

Eigen::MatrixXd transformSelenographicToCraterFrame(const Eigen::Vector3d& position, const Eigen::Matrix3d& T_e2m) {
  Eigen::Matrix3d h_m(3,3);
  Eigen::Vector3d u_north_pole = getNorthPoleUnitVector();
  // form transformation matrix
  // eq. 34 from Christian, Derksen, Watkins [2020]
  h_m << T_e2m.col(0), T_e2m.col(1), position;
  // express matrix in homogeneous form (eq. 40)
  Eigen::MatrixXd h_k(4,3);
  h_k << h_m.row(0), h_m.row(1), h_m.row(2), u_north_pole.transpose();
  return h_k;
}
