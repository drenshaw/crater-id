#include "quadrics.h"
#include "io.h"
#include "math_utils.h"
#include "conics.h"

#include <iomanip>

// This is the base constructor
// Technically, these are a form of quadrics (quadratic surfaces) known
// as a quadric disk, which is a flat circular disk, which we use to
// represent crater rims.
// TODO: Should I rename this to `QuadricDisk`, then?
Quadric::Quadric( const Eigen::Vector3d& position, 
                  const double radius, 
                  const Eigen::Vector3d& surface_normal, 
                  const std::string id) : 
                    id_(id),
                    radius_(radius),
                    surface_point_{position},
                    plane_(surface_normal.normalized(), position.norm()) {
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
  Eigen::Matrix3d T_e2m_ = getQuadricTransformationMatrix();
  // plane_ = SurfacePointToPlane(T_e2m_, surface_point_);
  plane_ = SurfacePointToPlane(T_e2m_, surface_point_);
}
Quadric::Quadric(const Eigen::Vector3d& pt1, const Eigen::Vector3d& pt2, const Eigen::Vector3d& pt3, std::string id) :
  id_(id),
  plane_(Eigen::Hyperplane<double, 3>::Through(pt1, pt2, pt3)) {
  // This is the 3d quadric disk formed by three 3D points
  Eigen::Hyperplane<double, 3> plane_ = Eigen::Hyperplane<double, 3>::Through(pt1, pt2, pt3);
  // Using ideas from github.com/sergarrido/random/blob/master/circle3d
  // math.stackexchange.com/questions/1076177/3d-coordinates-of-circle-center-given-three-point-on-the-circle
  
  Eigen::Vector3d p12 = pt2 - pt1;
  Eigen::Vector3d p13 = pt3 - pt1;
  double v11, v22, v12;
  v11 = p12.dot(p12);
  v22 = p13.dot(p13);
  v12 = p12.dot(p13);

  double base = 0.5/(v11*v22-v12*v12);
  double k1 = base*v22*(v11-v12);
  double k2 = base*v11*(v22-v12);
  Eigen::Vector3d center = pt1 + p12*k1 + p13*k2; // center
  
  // // Alternative solution
  // Eigen::Vector3d u = (B - A).normalized();
  // Eigen::Vector3d w = (C - A).cross(u).normalized();
  // Eigen::Vector3d v = w.cross(u);
  // double bx = (B - A).dot(u);
  // // double by = 0;
  // double cx = (C - A).dot(u);
  // double cy = (C - A).dot(v);
  // double num = std::pow((cx-bx/2), 2) + std::pow(cy, 2) - std::pow(bx/2, 2);
  // double h = num / (2*cy);
  // Eigen::Vector3d center = A + (bx/2)*u + h*v;

  this->surface_point_ = center;
  this->radius_ = (pt1 - center).norm();
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
        
Eigen::Matrix4d Quadric::generateQuadricLocus() const {
  return GenerateQuadricLocusFromRadiusNormal(surface_point_, radius_);
}

Eigen::Matrix3d Quadric::getQuadricTransformationMatrix() const {
  return getENUFrame(plane_.normal());
}

Eigen::Matrix4d Quadric::getLocus() const {
  return generateQuadricLocus();
}

Eigen::Vector3d Quadric::getLocation() const {
  return surface_point_;
}

void Quadric::getLocation(Eigen::Vector3d& location) const {
  location = surface_point_;
}

Eigen::Vector3d Quadric::getNormal() const {
  return plane_.normal();
}

void Quadric::getNormal(Eigen::Vector3d& surface_normal) const {
  surface_normal = plane_.normal();
}

Eigen::Hyperplane<double, 3> Quadric::getPlane() const {
  return plane_;
}

void Quadric::getPlane(Eigen::Hyperplane<double, 3>& hyperplane) const {
  hyperplane = plane_;
}

double Quadric::getAngleBetweenQuadrics(const Quadric& other_quadric) const {
  return getAngleBetweenVectors(this->getNormal(), other_quadric.getNormal());
}

Eigen::Vector3d Quadric::getAxisNormalToQuadrics(const Quadric& other_quadric) const {
  Eigen::Vector3d axis_ret(3);
  try {
    axis_ret = getAxisNormalToVectors(
      this->getNormal(), 
      other_quadric.getNormal());
  }
  catch (const std::exception& e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    throw std::runtime_error("Quadric normals are too close to form an axis of rotation.");
  }
  return axis_ret;
}

double Quadric::getRadius() const {
  return radius_;
}


std::ostream& operator<<(std::ostream& os, const Quadric& quad) {
  std::streamsize ss = std::cout.precision();
  std::streamsize sw = std::cout.width();
  return os 
    << "Quadric -> \tID: '" << quad.id_ << "'"
    << std::fixed << std::setw(8) << std::setprecision(1)
    << "\n\t\tLocation: (" << quad.surface_point_.transpose() << ") km | "
    << "\tRadius: " << quad.radius_ << "km"
    << "\n\t\tPlane: [" << quad.plane_.coeffs().transpose() << " ] "
    << std::setprecision(ss) << std::setw(sw) << std::defaultfloat
    ;
}

double calculateCraterRimFromRadius(const double radius) {
  return sqrt(pow(R_MOON, 2) - pow(radius, 2));
}

Eigen::Hyperplane<double, 3> SurfacePointToPlane(const Eigen::Matrix3d& T_e2m, 
                                    const Eigen::Vector3d& surface_point) {
  Eigen::Vector3d u_north_pole = getNorthPoleUnitVector();
  Eigen::Vector3d plane_normal = T_e2m * u_north_pole;
  double rho = surface_point.dot(plane_normal);
  Eigen::Hyperplane<double, 3> plane(plane_normal, -rho);
  return plane;
}

Eigen::Matrix4d GenerateQuadricLocusFromRadiusNormal(const Eigen::Vector3d& position, const double radius) {
  Conic conic(radius, radius, 0, 0, 0);
  Eigen::Matrix3d conic_envelope = conic.getEnvelope();

  Eigen::Matrix3d T_enu_to_ref = getENUFrame(position);
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
