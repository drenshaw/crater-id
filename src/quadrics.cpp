#include "quadrics.h"
#include "io.h"
#include "vector_math.h"
#include "conics.h"

Quadric::Quadric(const Eigen::Vector3d& position, const double radius, const Eigen::Vector3d& surface_normal, const std::string id)
{
    MakeQuadric(position, radius, surface_normal, id);
}
Quadric::Quadric(const Eigen::Vector3d& position, const double radius, const std::string id)
{
    MakeQuadric(position, radius, id);
}
Quadric::Quadric(const std::string id, const Eigen::Vector3d& position, const double radius) {
    MakeQuadric(position, radius, id);
}
Quadric::Quadric(const double lat, const double lon, const double radius, const std::string id)
{
    MakeQuadric(lat, lon, radius, id);
}
Quadric::Quadric(const std::string id, const double lat, const double lon, const double radius) {
    MakeQuadric(lat, lon, radius, id);
}
Quadric::Quadric(const std::string id, const Eigen::Vector3d& position, const double radius, const Eigen::Vector3d& surface_normal) {
    MakeQuadric(position, radius, surface_normal, id);
}
void Quadric::MakeQuadric(const Eigen::Vector3d& position, const double radius, const Eigen::Vector3d& surface_normal, const std::string id) {
    surface_point_ = position;
    assert(("The surface normal is not defined or ==0.", surface_normal.norm() != 0));
    surface_normal_ = surface_normal;
    surface_normal_.normalize();
    radius_ = radius;
    id_ = id;
    T_e2m_ = GetQuadricTransformationMatrix();
    std::cout << "T_e2m: (determinant " << T_e2m_.determinant() << ")\n" << T_e2m_ << std::endl;
    Conic conic(radius, radius, 0, 0, 0);
    locus_ = GenerateQuadricLocus();
    std::cout << "Locus \n" << locus_ << std::endl;
    envelope_ = getAdjugateMatrix(locus_);
    std::tie(plane_, plane_normal_) = SurfacePointToPlane(T_e2m_, surface_point_);
}
// assumes that the surface normal of the crater is in line with the position wrt Moon center
void Quadric::MakeQuadric(const Eigen::Vector3d& position, const double radius, const std::string id) {
    Eigen::Vector3d surface_normal = position.normalized();
    MakeQuadric(position, radius, surface_normal, id);
}
void Quadric::MakeQuadric(const double lat, const double lon, const double radius, const std::string id) {
    if(radius > R_MOON) {
        std::cerr << "Invalid crater radius " << radius << "exceeds body radius\n";
    }
    double r_crater_bodycenter = calculateCraterRimFromRadius(radius);
    Eigen::Vector3d position = r_crater_bodycenter*latlon2bearing(lat, lon);
    MakeQuadric(position, radius, id);
}
        
Eigen::Matrix4d Quadric::GenerateQuadricLocus() {
    Eigen::Matrix4d quadric_locus;
    return GenerateQuadricFromRadiusNormal(surface_point_, radius_);
}

Eigen::Matrix3d Quadric::GetQuadricTransformationMatrix() {
    return getENUFrame(surface_normal_);
}

Eigen::Matrix4d Quadric::GetLocus() {
    return locus_;
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

std::tuple<Eigen::Vector4d, Eigen::Vector3d> SurfacePointToPlane(const Eigen::Matrix3d& T_e2m, 
                                                                 const Eigen::Vector3d& surface_point) {
    Eigen::Vector3d u_north_pole = getNorthPoleUnitVector();
    Eigen::Vector3d plane_normal = T_e2m * u_north_pole;
    double rho = surface_point.dot(plane_normal);
    Eigen::Vector4d plane;
    plane << plane_normal, -rho;
    return {plane, plane_normal};
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
