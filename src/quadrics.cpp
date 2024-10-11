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
                    center_{position},
                    plane_(surface_normal.normalized(), position.norm()
                  ) {
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
  plane_ = SurfacePointToPlane(T_e2m_, center_);
}
Quadric::Quadric(const Eigen::Vector3d& pt1, const Eigen::Vector3d& pt2, const Eigen::Vector3d& pt3, std::string id) :
  id_(id),
  plane_(Eigen::Hyperplane<double, 3>::Through(pt1, pt2, pt3)) {
  // This is the 3d quadric disk formed by three 3D points
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

  this->center_ = center;
  this->radius_ = (pt1 - center).norm();
}

// These are delegating constructors
Quadric::Quadric(const Eigen::Vector3d& position, const double radius, const std::string id) : 
  Quadric(position, radius, position.normalized(), id) {}
Quadric::Quadric(const std::string id, const Eigen::Vector3d& position, const double radius) :
  Quadric(position, radius, id) {}
Quadric::Quadric(const double lat, const double lon, const double radius, const std::string id) :
  Quadric(latlonrad2XYZ(lat, lon, radius), radius, id) {}
Quadric::Quadric(const std::string id, const double lat, const double lon, const double radius) :
  Quadric(latlonrad2XYZ(lat, lon, radius), radius, id) {}
Quadric::Quadric(const std::string id, const Eigen::Vector3d& position, const double radius, const Eigen::Vector3d& surface_normal) :
  Quadric(position, radius, surface_normal, id) {}

bool Quadric::operator==(const Quadric& other_quadric) const {
  return  // TODO: Sometimes the normal can be opposite direction; check
    isSamePlane(this->getPlane(),other_quadric.getPlane()) &&
    this->getLocation().isApprox(other_quadric.getLocation()) &&
    std::abs(this->getRadius() - other_quadric.getRadius()) < 1e-3;
}

bool Quadric::operator!=(const Quadric& other_quadric) const {
  return  !operator==(other_quadric);
}

bool Quadric::operator==(const Quadric* other_quadric) const {
  return  operator==(*other_quadric);
}

bool Quadric::operator!=(const Quadric* other_quadric) const {
  return  !operator==(other_quadric);
}

Eigen::Matrix4d Quadric::generateQuadricLocusFromPointRadius() const {
  return GenerateQuadricLocus(this->center_, this->radius_);
}

Eigen::Matrix4d Quadric::generateQuadricEnvelopeFromPointRadius() const {
  return GenerateQuadricEnvelope(this->center_, this->radius_);
}

Eigen::Matrix3d Quadric::getQuadricTransformationMatrix() const {
  return getENUFrame(this->plane_.normal());
}

Eigen::Matrix4d Quadric::getLocus() const {
  return this->generateQuadricLocusFromPointRadius();
}

Eigen::Matrix4d Quadric::getEnvelope() const {
  return this->generateQuadricEnvelopeFromPointRadius();
}

Eigen::Vector3d Quadric::getLocation() const {
  return this->center_;
}

void Quadric::getLocation(Eigen::Vector3d& location) const {
  location = this->center_;
}

Eigen::Vector3d Quadric::getNormal() const {
  return this->plane_.normal();
}

void Quadric::getNormal(Eigen::Vector3d& surface_normal) const {
  surface_normal = this->plane_.normal();
}

Eigen::Hyperplane<double, 3> Quadric::getPlane() const {
  return this->plane_;
}

void Quadric::getPlane(Eigen::Hyperplane<double, 3>& hyperplane) const {
  hyperplane = this->plane_;
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
    throw std::runtime_error("Quadric normals are nearly parallel.");
  }
  return axis_ret;
}

double Quadric::getRadius() const {
  return this->radius_;
}

std::string Quadric::getID() const {
  return this->id_;
}

Conic Quadric::projectToConic(const Eigen::MatrixXd& proj_mtx) const{
  Eigen::Matrix3d conic_envelope = proj_mtx * this->getEnvelope() * proj_mtx.transpose();
  Conic conic(adjugate(conic_envelope));
  return conic;
}


std::ostream& operator<<(std::ostream& os, const Quadric& quad) {
  std::streamsize ss = std::cout.precision();
  std::streamsize sw = std::cout.width();
  return os 
    << "Quadric -> \tID: '" << quad.id_ << "'"
    << std::fixed << std::setw(8) << std::setprecision(2)
    << "\n\t\tLocation: (" << quad.center_.transpose() << ") km | "
    << "\tRadius: " << quad.radius_ << "km"
    << "\n\t\tPlane: [" << quad.plane_.coeffs().transpose() << " ] "
    << std::setprecision(ss) << std::setw(sw) << std::defaultfloat
    ;
}


/***************UTILS**************/

bool isSamePlane(const Eigen::Hyperplane<double, 3>& p1, const Eigen::Hyperplane<double, 3>& p2, const double thresh) {
  if(p1.isApprox(p2))
    return true;
  if((p1.normal() - p2.normal()).norm() < 1e-3 && std::abs(p1.offset() - p2.offset()) < thresh)
    return true;
  if((p1.normal() + p2.normal()).norm() < 1e-3 && std::abs(p1.offset() + p2.offset()) < thresh)
    return true;
  return false;
}

bool isSamePlane(const Quadric& quad1, const Quadric& quad2, const double thresh) {
  return isSamePlane(quad1.getPlane(), quad2.getPlane(), thresh);
}

Eigen::Vector3d latlonrad2XYZ(const double lat, const double lon, const double radius) {
  const double dist = calculateCraterRimFromRadius(radius);
  return dist * latlon2bearing(lat, lon);
}

double calculateCraterRimFromRadius(const double radius) {
  return std::sqrt(std::pow(R_MOON, 2) - std::pow(radius, 2));
}

Eigen::Hyperplane<double, 3> SurfacePointToPlane(const Eigen::Matrix3d& T_e2m, 
                                    const Eigen::Vector3d& surface_point) {
  Eigen::Vector3d u_north_pole = getNorthPoleUnitVector();
  Eigen::Vector3d plane_normal = T_e2m * u_north_pole;
  double rho = surface_point.dot(plane_normal);
  Eigen::Hyperplane<double, 3> plane(plane_normal, -rho);
  return plane;
}

Eigen::Matrix4d GenerateQuadricLocus(const Eigen::Vector3d& position, const double radius) {
  Eigen::Matrix4d quadric_envelope = GenerateQuadricEnvelope(position, radius);
  Eigen::Matrix4d quadric_locus = adjugate(quadric_envelope);
  // std::cout << __func__ << ": Locus:\n" << quadric_locus/quadric_locus(3,3) << std::endl;
  // bool success = normalizeDeterminant(quadric_locus);
  return quadric_locus;
}

Eigen::Matrix4d GenerateQuadricEnvelope(const Eigen::Vector3d& position, const double radius) {
  Conic conic(radius, radius, 0, 0, 0);
  Eigen::Matrix3d conic_envelope = conic.getEnvelope();
  Eigen::Matrix3d T_enu_to_ref = getENUFrame(position);
  Eigen::MatrixXd h_k = transformSelenographicToCraterFrame(position, T_enu_to_ref);
  // std::cout << __func__ << ": h_k:\n" << h_k << std::endl;
  // eq 40 of Christian, Derksen, and Watkins [2020]
  Eigen::Matrix4d quadric_envelope = ConicEnvelopeToQuadricEnvelope(conic_envelope, h_k);
  // std::cout << __func__ << ": Envelope:\n" << quadric_envelope/quadric_envelope(3,3) << std::endl;
  return quadric_envelope;
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
  h_k << h_m, u_north_pole.transpose();
  return h_k;
}

Eigen::Matrix4d makeSphere(const double radius) {
  double recip_radius = 1/std::pow(radius,2.0);
  Eigen::Matrix4d sphere = recip_radius * Eigen::Matrix4d::Identity();
  sphere(3,3) = -1.0;
  return sphere;
}

Eigen::Matrix4d makeEllipsoid(const Eigen::Vector3d& radii) {
  Eigen::Vector3d recip_radii = radii.array().pow(-2);
  Eigen::Vector4d diag;
  diag << recip_radii, -1;
  Eigen::Matrix4d sphere = Eigen::Matrix4d::Identity();
  sphere.diagonal() << diag;
  return sphere;
}
