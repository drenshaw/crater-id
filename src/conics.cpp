#include "conics.h"
#include "math_utils.h"

#include <iostream>
// #include <stdexcept>
#include <tuple>
#include <array>
#include <cmath>
#include <random>
#include <eigen3/Eigen/Dense>
#include <optional>

int Conic::next_id = 0;

// Conic::Conic() : Conic(1, 1, 0, 0, 0) {}
Conic::Conic(const double semimajor_axis, 
             const double semiminor_axis, 
             const double x_center, 
             const double y_center, 
             const double angle) : 
    semimajor_axis_{semimajor_axis},
    semiminor_axis_{semiminor_axis},
    x_center_      {x_center},
    y_center_      {y_center},
    angle_         {angle} {
  setID();
}
Conic::Conic(const std::array<double, GEOMETRIC_PARAM>& geom_arr) : 
    semimajor_axis_{geom_arr.at(0)},
    semiminor_axis_{geom_arr.at(1)},
    x_center_      {geom_arr.at(2)},
    y_center_      {geom_arr.at(3)},
    angle_         {geom_arr.at(4)} {
  setID();
}
Conic::Conic(const std::vector<double>& geom_vec) : 
    semimajor_axis_{geom_vec.at(0)},
    semiminor_axis_{geom_vec.at(1)},
    x_center_      {geom_vec.at(2)},
    y_center_      {geom_vec.at(3)},
    angle_         {geom_vec.at(4)} {
  this->setID();
}
Conic::Conic(const std::array<double, IMPLICIT_PARAM>& impl) {
  setImplicitParameters(impl);
}
Conic::Conic(const Eigen::Matrix3d& locus) {
  this->setLocus(locus);
}
Conic::Conic(const std::vector<Eigen::Vector2d>& points) {
  if(points.size() < 5) {
    std::cerr << "Too few points to form an ellipse: " 
              << "need >=5, given " << points.size() << std::endl;
    throw std::runtime_error("Too few points to form an ellipse.");              
  }
  std::array<double, IMPLICIT_PARAM> impl;
  impl = ellipseFitLstSq(points);
  setImplicitParameters(impl);
}

void Conic::setID() {
  this->id_ = next_id++;
}

int Conic::getID() const {
  return this->id_;
}

bool Conic::operator==(const Conic& other_conic) const {
  return  almost_equal(this->semimajor_axis_, other_conic.getSemiMajorAxis()) &&
          almost_equal(this->semiminor_axis_, other_conic.getSemiMinorAxis()) &&
          almost_equal(this->x_center_, other_conic.getCenterX()) &&
          almost_equal(this->y_center_, other_conic.getCenterY()) &&
          almost_equal(this->angle_, other_conic.getAngle());
}

bool Conic::operator!=(const Conic& other_conic) const {
  return  !operator==(other_conic);
}

bool Conic::operator==(const Conic* other_conic) const {
  return  operator==(*other_conic);
}

bool Conic::operator!=(const Conic* other_conic) const {
  return  !operator==(other_conic);
}

std::ostream& operator<<(std::ostream& os, const Conic& conic) {
    return os 
        << "Conic-> \tID: '" << conic.id_ << "'"
        << "\n\tSemi axes: " 
        << conic.getSemiMajorAxis() << " | "
        << conic.getSemiMinorAxis()
        << "\n\tCenter: (" << conic.getCenterX() << " , "
        << conic.getCenterY() << ") "
        << "\n\tAngle (deg): " << conic.getAngleDeg() << " ";
}

void Conic::setGeometricParameters(const std::array<double, GEOMETRIC_PARAM>& geom_vec) {
  setGeometricParameters( geom_vec.at(0), 
                          geom_vec.at(1), 
                          geom_vec.at(2), 
                          geom_vec.at(3), 
                          geom_vec.at(4));
}

void Conic::setGeometricParameters(const std::vector<double>& geom_vec) {
  setGeometricParameters( geom_vec.at(0), 
                          geom_vec.at(1), 
                          geom_vec.at(2), 
                          geom_vec.at(3), 
                          geom_vec.at(4));
}

void Conic::setGeometricParameters( const double semimajor_axis, 
                                    const double semiminor_axis, 
                                    const double x_center, 
                                    const double y_center, 
                                    const double angle) {
  if(std::isnan(semimajor_axis) || std::isnan(semiminor_axis) || 
     std::isnan(x_center)       || std::isnan(y_center) || 
     std::isnan(angle)) {
    throw std::runtime_error("One or more parameters is NaN");
  }
  semimajor_axis_ = semimajor_axis;
  semiminor_axis_ = semiminor_axis;
  x_center_ = x_center;
  y_center_ = y_center;
  angle_ = wrap_npi2_pi2(angle);
}

void Conic::setSemimajorAxis(const double semimajor_axis) {
  semimajor_axis_ = semimajor_axis;
}

void Conic::setSemiminorAxis(const double semiminor_axis) {
  semiminor_axis_ = semiminor_axis;
}

void Conic::setCenterX(const double x_center) {
  x_center_ = x_center;
}

void Conic::setCenterY(const double y_center) {
  y_center_ = y_center;
}

void Conic::setAngle(const double angle) {
  angle_ = angle;
}


void Conic::setImplicitParameters(const std::array<double, IMPLICIT_PARAM>& impl_params) {
  std::optional<std::array<double, GEOMETRIC_PARAM> > geom_params = implicit2Geom(impl_params);
  if(geom_params.has_value()) {
    setGeometricParameters(geom_params.value());
  }
  else {
    throw std::runtime_error("Invalid implicit parameters");
  }
}

void Conic::setLocus(const Eigen::Matrix3d& locus) {
  std::optional<std::array<double, GEOMETRIC_PARAM> > geom_params;
  geom_params = fromLocus(locus);
  if(geom_params.has_value()) {
    setGeometricParameters(geom_params.value());
  }
  else {
    throw std::runtime_error("Provided locus is invalid");
  }
}

std::array<double, GEOMETRIC_PARAM> Conic::getGeom() const {
  return {semimajor_axis_, semiminor_axis_, x_center_, y_center_, angle_};
}

std::array<double, IMPLICIT_PARAM> Conic::getImplicit() const {
  return toImplicit();
}

Eigen::Matrix3d Conic::getLocus() const {
  return toLocus();
}

Eigen::Matrix3d Conic::getEnvelope() const {
  return adjugate(toLocus());
}

Eigen::Matrix3d Conic::toLocus() const {
  return implicit2Locus(this->getImplicit());
}

std::optional<std::array<double, GEOMETRIC_PARAM> > Conic::fromLocus(const Eigen::Matrix3d& locus) const {
  return locus2Geom(locus);
}


std::array<double, IMPLICIT_PARAM> Conic::toImplicit() const {
  return geom2Implicit( this->semimajor_axis_, 
                        this->semiminor_axis_, 
                        this->x_center_, 
                        this->y_center_, 
                        this->angle_);
}

bool Conic::intersectsConic(const Eigen::Matrix3d& Aj, 
                            std::tuple<Eigen::Vector3d, Eigen::Vector3d>& gh) const {
  Eigen::Matrix3d Ai = this->getLocus();
  return invariants::intersectionLines(Ai, Aj, gh);
}

bool Conic::intersectsConicLines(const Conic& conicB, 
                                   std::tuple<Eigen::Vector3d, Eigen::Vector3d>& gh) const {
  Eigen::Matrix3d Ai = this->getLocus();
  Eigen::Matrix3d Aj = conicB.getLocus();
  return invariants::intersectionLines(Ai, Aj, gh);
}

double Conic::getCenterX() const {
  return x_center_;
}

double Conic::getCenterY() const {
  return y_center_;
}

Eigen::Vector2d Conic::getCenter() const {
  Eigen::Vector2d center;
  center << x_center_, y_center_;
  return center;
}

void Conic::getCenter(cv::Point2d& center) const {
  center.x = x_center_;
  center.y = y_center_;
}

void Conic::getCenter(Eigen::Vector2d& center) const {
  center << x_center_, y_center_;
}

double Conic::getSemiMajorAxis() const {
  return semimajor_axis_;
}

double Conic::getSemiMinorAxis() const {
  return semiminor_axis_;
}

cv::Size Conic::getSemiAxes() const {
  cv::Size axes(semimajor_axis_, semiminor_axis_);
  return axes;
}

void Conic::getSemiAxes(Eigen::Vector2d& semiaxes) const {
  semiaxes << semimajor_axis_, semiminor_axis_;
}

void Conic::getSemiAxes(cv::Point& semiaxes) const {
  semiaxes.x = semimajor_axis_;
  semiaxes.y = semiminor_axis_;
}

void Conic::getSemiAxes(cv::Size& semiaxes) const {
  semiaxes.width  = semimajor_axis_;
  semiaxes.height = semiminor_axis_;
}

// GetSize is an alias for GetSemiAxes
cv::Size Conic::getSize() const {
  return getSemiAxes();
}

void Conic::getSize(cv::Size& semiaxes) const {
  return getSemiAxes(semiaxes);
}

double Conic::getAngle() const {
  return angle_;
}

double Conic::getAngleDeg() const {
  return rad2deg(angle_);
}

void Conic::normalizeImplicitParams() {
  // normalizeImplicitParameters(this->impl)
}

bool Conic::chooseIntersection(const Conic& other, Eigen::Vector3d& l) const {
  std::tuple<Eigen::Vector3d, Eigen::Vector3d> gh;
  if(!invariants::intersectionLines(getLocus(), other.getLocus(), gh)) {
    std::cerr<<"IntersectionLines error ij\n";
    return false;
  }
  if(!invariants::chooseIntersection(gh, getCenter(), other.getCenter(), l)) {
    std::cerr<<"ChooseIntersection error jk\n";
    return false;
  }

  return true;
}

double Conic::getEccentricity() const {
  double a = this->getSemiMajorAxis();
  double b = this->getSemiMinorAxis();
  return std::sqrt(1 - std::pow(b, 2)/std::pow(a, 2));
}


/*********************************************************/
/***********************Conic Utils***********************/
/*********************************************************/

std::optional<std::array<double, IMPLICIT_PARAM> > locus2Implicit(const Eigen::Matrix3d& locus) {
  Eigen::Matrix3d loc = locus;
  std::optional<std::array<double, IMPLICIT_PARAM> > params;
  normalizeDeterminant(loc);
  // loc = locus(2,2) == 0 ? locus : locus/locus(2,2);
  
  if(loc.determinant() == 0) {
    std::cerr << __func__ << "--> Locus is singular:\n" << loc << std::endl;
    return params;
  }
  /*
      |  a   b/2  d/2 |
  C = | b/2   c   e/2 |
      | d/2  e/2   f  |
  */
  if(loc.hasNaN()) {
    std::cerr << __func__ << "--> Locus has NaN:\n" << loc << std::endl;
    throw std::runtime_error("The above locus has NaN value");
  }
  if(!isEllipse(loc)) {
    std::cerr << __func__ << "--> \n" << loc << std::endl;
    throw std::runtime_error("The above input matrix does not represent an ellipse.");
  }
  double A, B, C, D, E, F;
  A =    loc.coeff(0,0);
  B =  2*loc.coeff(0,1);
  C =    loc.coeff(1,1);
  D =  2*loc.coeff(0,2);
  E =  2*loc.coeff(1,2);
  F =    loc.coeff(2,2); 
  params = {A, B, C, D, E, F};
  return params;
}

std::optional<std::array<double, GEOMETRIC_PARAM> > locus2Geom(const Eigen::Matrix3d& locus) {
  const std::optional<std::array<double, IMPLICIT_PARAM> > impl_params = locus2Implicit(locus);
  if(impl_params.has_value())
    return implicit2Geom(impl_params.value());
  else {
    return std::nullopt;
  }
}

bool isEllipse(const Eigen::Matrix3d& conic_locus) {
    // Check if the determinant of the top-left 2x2 submatrix is positive
    Eigen::Matrix2d subMatrix = conic_locus.topLeftCorner<2, 2>();
    return subMatrix.determinant() > 0;
}

// void extractEllipseParameters(const Eigen::Matrix3d& A, double& semiMajor, double& semiMinor, Eigen::Vector2d& center, double& angle) {
//     // Extract the top-left 2x2 submatrix
//     Eigen::Matrix2d subMatrix = A.topLeftCorner<2, 2>();
    
//     // Eigenvalue decomposition
//     Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(subMatrix);
//     if (eigensolver.info() != Eigen::Success) {
//         std::cerr << "Eigenvalue decomposition failed!" << std::endl;
//         return;
//     }
    
//     // Eigenvalues and eigenvectors
//     Eigen::Vector2d eigenvalues = eigensolver.eigenvalues();
//     Eigen::Matrix2d eigenvectors = eigensolver.eigenvectors();
    
//     // Semi-major and semi-minor axes
//     semiMajor = 1.0 / std::sqrt(eigenvalues.minCoeff());
//     semiMinor = 1.0 / std::sqrt(eigenvalues.maxCoeff());
    
//     // Center of the ellipse
//     center = -subMatrix.inverse() * A.topRightCorner<2, 1>();
    
//     // Angle with respect to the x-axis
//     angle = std::atan2(eigenvectors(1, 0), eigenvectors(0, 0));
// }

std::optional<std::array<double, GEOMETRIC_PARAM> > implicit2Geom(const std::array<double, IMPLICIT_PARAM>& impl_params) {
  // double numerator, denominator_a, denominator_b;
  std::optional<std::array<double, GEOMETRIC_PARAM> > geom;
  double B2_minus_4AC, xc, yc, phi;
  double semimajor_axis, semiminor_axis;//, amc2, b2;
  /*
      |  A   B/2  D/2 |
  C = | B/2   C   E/2 |
      | D/2  E/2   F  |
  */
  if(std::any_of(impl_params.begin(), impl_params.end(), 
                                   [](double x) { return std::isnan(x); })) {
    std::cout << __func__ << " -- NaN value found\n\t";
    for(const auto& elem : impl_params) {
      std::cout << elem << ", ";
    }
    std::cout << std::endl;
    return geom;
  }
  double A = impl_params.at(0);
  double B = impl_params.at(1);
  double C = impl_params.at(2);
  double D = impl_params.at(3);
  double E = impl_params.at(4);
  double F = impl_params.at(5);
  // Compute discriminant for quadratic equation
  B2_minus_4AC = std::pow(B, 2) - 4*A*C;
  assert(B2_minus_4AC < 0);
  // Compute ellipse center (See Eq 4.16 in [Christian, 2010])
  xc = (2*C*D - B*E) / B2_minus_4AC;
  yc = (2*A*E - B*D) / B2_minus_4AC;
  // Compute ellipse semimajor axis (a) and seminor axis (b)
  // (See Eqs 4.17 and 4.18 in [Christian, 2010])
  double num;
  num = 2*((A*E*E + C*D*D - B*D*E)/(4*A*C - B*B) - F);
  double den = std::sqrt(std::pow(A-C,2) + std::pow(B,2));
  assert(!std::isnan(den));
  // assert(den != 0 && !std::isnan(den));
  // TODO: I think this departs from Christian 2010; find source
  semimajor_axis = std::sqrt(num/(A+C-den));
  semiminor_axis = std::sqrt(num/(A+C+den));

  // Compute angle from the x axis to semimajor axis direction
  // (See Eq 4.19 in [Christian, 2010])
  // NOTE: do NOT use `atan2` for phi calculation
  phi = (B == 0) ? 0.0 : 0.5*std::atan(B / (A - C));
  if(A > C)
    phi += M_PI_2;
  if(semimajor_axis < semiminor_axis) {
    // std::cout << "Swappin' axes\n";
    std::swap(semimajor_axis, semiminor_axis);
    assert(semimajor_axis>=semiminor_axis);
    phi += M_PI_2;
  }
  // TODO: figure out why the negative angle is visually correct
  phi = wrap_npi2_pi2(phi);
  // auto roundd = [](double val) {return std::round(val*1e3)/1e3;};
  // auto roundd = [](double val) {return val;};
  if(std::isnan(semimajor_axis) || std::isnan(semiminor_axis)) {
    throw std::runtime_error("At least one axis is NaN");
  }
  geom = {
    semimajor_axis, 
    semiminor_axis, 
    xc, 
    yc, 
    phi};
  if(vectorContainsNaN(geom.value())) {
    std::cerr 
      << __func__ << "-->\n"
      << "Semimajor axis : " << semimajor_axis << std::endl
      << "Semiminor axis : " << semiminor_axis << std::endl
      << "Center (x)     : " << xc << std::endl
      << "Center (y)     : " << yc << std::endl
      << "Angle (phi)    : " << phi << std::endl;
    throw std::runtime_error("Above geometric parameters contain NaN value(s).");
  }
  return geom;
}

std::array<double, IMPLICIT_PARAM> geom2Implicit( const double semimajor_axis, 
                                                  const double semiminor_axis, 
                                                  const double x_center, 
                                                  const double y_center, 
                                                  const double angle) {
  // unpack geometric parameters
  double a   = semimajor_axis;
  double b   = semiminor_axis;
  double xc  = x_center;
  double yc  = y_center;
  double phi = angle;

  // perform some computations beforehand
  double a2 = pow(a, 2);
  double b2 = pow(b, 2);
  double sin_phi = sin(phi);
  double cos_phi = cos(phi);
  double xc_2 = pow(xc, 2);
  double yc_2 = pow(yc, 2);
  double sin_phi_2 = pow(sin_phi, 2);
  double cos_phi_2 = pow(cos_phi, 2);

  std::array<double, IMPLICIT_PARAM> coeff;

  // Using conversion found in Christion 2020
  // "Lunar crater identification in digital images"
  coeff.at(0) =  a2*sin_phi_2 + b2*cos_phi_2;
  coeff.at(1) =  2*(b2-a2)*sin_phi*cos_phi;
  coeff.at(2) =  a2*cos_phi_2 + b2*sin_phi_2;
  coeff.at(3) = -2*coeff[0]*xc - coeff[1]*yc;
  coeff.at(4) = -coeff[1]*xc - 2*coeff[2]*yc;
  coeff.at(5) =  coeff[0]*xc_2 + coeff[1]*xc*yc + coeff[2]*yc_2 - a2*b2;

  normalizeImplicitParameters(coeff);
  return coeff;
}

Eigen::Matrix3d implicit2Locus(const std::array<double, IMPLICIT_PARAM>& impl_params) {
  Eigen::Matrix3d locus_mtx(3,3);
  double A, B, C, D, E, F;

  F = impl_params.at(5);
  double F_recip = 1/F;
  A = F_recip * impl_params.at(0);
  B = F_recip * impl_params.at(1);
  C = F_recip * impl_params.at(2);
  D = F_recip * impl_params.at(3);
  E = F_recip * impl_params.at(4);
  locus_mtx << A,  B/2, D/2,
              B/2,  C,  E/2,
              D/2, E/2, 1.0;

  return locus_mtx;
}

void normalizeImplicitParameters(std::array<double, IMPLICIT_PARAM>& impl_params) {
  // double vecNormRecip = 1/vectorNorm(impl_params);
  double vecNormRecip = 1/impl_params.at(5);
  std::transform(impl_params.begin(), impl_params.end(), impl_params.begin(),
               std::bind(std::multiplies<double>(), std::placeholders::_1, vecNormRecip));
}

// TODO: possibly create template or call to same function for std::array above
void normalizeImplicitParameters(std::vector<double>& impl_params) {
  // double vecNormRecip = 1/vectorNorm(impl_params);
  double vecNormRecip = 1/impl_params.at(5);
  std::transform(impl_params.begin(), impl_params.end(), impl_params.begin(),
               std::bind(std::multiplies<double>(), std::placeholders::_1, vecNormRecip));
}

std::array<double, IMPLICIT_PARAM> ellipseFitLstSq(const Eigen::MatrixXd& points) {
  std::array<double, IMPLICIT_PARAM> impl;
  int n_pts = points.rows();
  if(n_pts < 5) {
    impl.fill({std::nan("")});
  }
  Eigen::ArrayXd X = points.col(0).array();
  Eigen::ArrayXd Y = points.col(1).array();
  Eigen::MatrixXd Ax(n_pts, 5);
  Ax << X.pow(2), X*Y, Y.pow(2), X, Y;
  Eigen::MatrixXd bx = Eigen::MatrixXd::Ones(n_pts, 1);
  // Ax.template Eigen::BDCSVD<Eigen::ComputeThinU | Eigen::ComputeThinV>().solve(bx);
  Eigen::BDCSVD<Eigen::MatrixXd> SVD(Ax, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd solu = SVD.solve(bx);
  impl = {solu(0), solu(1), solu(2), solu(3), solu(4), -1.0};
  // Conic conic(impl);
  // // auto sol = Ax.colPivHouseholderQr().solve(bx);
  // Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(Ax.rows(), Ax.cols());
  // cod.setThreshold(Eigen::Default);
  // auto se = cod.compute(bx);
  return impl;
}

std::array<double, IMPLICIT_PARAM> ellipseFitLstSq(const std::vector<Eigen::Vector2d>& points) {
  std::array<double, IMPLICIT_PARAM> impl;
  int n_pts = points.size();
  if(n_pts < 5) {
    impl.fill({std::nan("")});
  }
  Eigen::MatrixXd pts_cam(n_pts, 2);
  std::vector<Eigen::Vector2d>::const_iterator it;
  for (it = points.begin(); it != points.end(); it++) {
    int index = std::distance(points.begin(), it);
    pts_cam(index, Eigen::all) = *it;
  }
  return ellipseFitLstSq(pts_cam);
}

void addNoise(const double mean, const double st_dev, std::vector<Eigen::Vector2d>& points) {
  if(mean == 0 && st_dev == 0) {
    return;
  }
  std::random_device rd;
  std::mt19937 gen(rd());
  // std::mt19937 gen(50);
  std::normal_distribution<double> dist(mean, st_dev); // Mean 0, standard deviation 1

  // Create an Eigen matrix and fill it with noise
  Eigen::MatrixXd A(3, 3);
  // for (int i = 0; i < A.rows(); ++i) {
  for(std::vector<Eigen::Vector2d>::iterator it = points.begin(); it != points.end(); it++) {
    // std::cout << "Points: " << (*it).transpose();
    (*it)(0) += dist(gen);
    (*it)(1) += dist(gen);
    // std::cout << "\tAdded noise: " << (*it).transpose() << std::endl;
  }
}

/*********************************************************/
/***********************INVARIANTS************************/
/*********************************************************/
namespace invariants {
Eigen::Matrix3d stackVectorsInto3x3Matrix(const Eigen::Vector3d& A, 
                                          const Eigen::Vector3d& B, 
                                          const Eigen::Vector3d& C) {
  Eigen::Matrix3d out_matrix;
  out_matrix.col(0) = A;
  out_matrix.col(1) = B;
  out_matrix.col(2) = C;
  return out_matrix;
}

double crossRatio(const Eigen::Vector3d& ref, 
                  const Eigen::Vector3d& A, 
                  const Eigen::Vector3d& B, 
                  const Eigen::Vector3d& C, 
                  const Eigen::Vector3d& D) {
  double numerator1, numerator2, denominator1, denominator2;
  numerator1 = stackVectorsInto3x3Matrix(A, B, ref).determinant();
  numerator2 = stackVectorsInto3x3Matrix(C, D, ref).determinant();
  denominator1 = stackVectorsInto3x3Matrix(A, C, ref).determinant();
  denominator2 = stackVectorsInto3x3Matrix(B, D, ref).determinant();
  assert(std::abs(denominator1)>1e-10);
  assert(std::abs(denominator2)>1e-10);
  double xratio = numerator1*numerator2/(denominator1*denominator2);
  return xratio;
}

bool IntersectConics(const Eigen::Matrix3d& Ai, 
                     const Eigen::Matrix3d& Aj, 
                     const double eig,
                     std::tuple<Eigen::Vector3d, Eigen::Vector3d>& gh) {
  double eps = 1e-16;
  Eigen::Matrix3d Bij(3, 3), Bij_star(3, 3);
  Eigen::Vector3d bkk_eig(3);
  Bij = eig*Ai + Aj;
  Bij_star = adjugate(Bij);

  bkk_eig = Bij_star.diagonal();
  auto valid = bkk_eig.array() < eps;
  if(!valid.all()) {
    return false;
  }
  Eigen::Vector3d::Index min_idx;
  min_idx = bkk_eig.minCoeff();

  // eq 87
  Eigen::Vector3d z = -Bij_star.row(min_idx)/sqrt(-bkk_eig[min_idx]);
  Eigen::Matrix3d D = Bij + crossMatrix(z);
  Eigen::Matrix3d Dabs = D.cwiseAbs();
  
  Eigen::Matrix3d::Index maxRow,maxCol;
  Dabs.maxCoeff(&maxRow, &maxCol);

  // pg 44 of Christian, Derksen, Watkins [2020]
  // g is non-zero column of D, h is non-zero row of D
  Eigen::Vector3d g = D.col(maxCol);
  Eigen::Vector3d h = D.row(maxRow);
  
  if ( vectorContainsNaN(g) ) {
      std::cerr << "g vector contains NaN's: " << g << std::endl;
      return false;
  }
  if ( vectorContainsNaN(h) ) {
      std::cerr << "h vector contains NaN's: " << h << std::endl;
      return false;
  }
  gh = {g, h};
  return true;
}

bool intersectionLines( const Eigen::Matrix3d& Ai, 
                        const Eigen::Matrix3d& Aj,
                        std::tuple<Eigen::Vector3d, Eigen::Vector3d>& gh) {
  // eq 77
  Eigen::Matrix3d combined = Aj*-Ai.inverse();
  // TODO: try suggested x = A.ldlt().solve(b) from section 2.12 of
  // https://geophydog.cool/post/eigen3_operations/ 
  Eigen::EigenSolver<Eigen::Matrix3d> eigensolver(combined);
  if (eigensolver.info() != Eigen::Success) {
    std::cerr << "Eigenvalue decomposition failed!" << std::endl;
    return false;
  }
  Eigen::Vector3cd eigenvalues = eigensolver.eigenvalues();
  // Here, using std::vector is appropriate since we don't know the length
  Eigen::Vector3d g, h;
  for(const std::complex<double>& eigen_elem : eigenvalues) {
    if(eigen_elem.imag() == 0 && IntersectConics(Ai, Aj, eigen_elem.real(), gh)) {
      std::tie(g, h) = gh; // TODO: What are we doing with these variables (g, h)?
      return true;
    }
  }
  return false;
}

bool chooseIntersection(const std::tuple<Eigen::Vector3d, Eigen::Vector3d>& gh, 
                        const Eigen::Vector2d& centerA, 
                        const Eigen::Vector2d& centerB,
                        Eigen::Vector3d& l) {
  Eigen::Vector3d g, h;
  std::tie(g, h) = gh;
  // convert centers to homogeneous coordinates
  Eigen::Vector3d centerAHom = centerA.homogeneous();
  Eigen::Vector3d centerBHom = centerB.homogeneous();
  // get line connecting the two centers
  Eigen::Vector3d lineOfCenters = centerAHom.cross(centerBHom);
  // get point where lineOfCenters and g intersect
  Eigen::Vector2d gIntersect = g.cross(lineOfCenters).hnormalized();
  // get point where lineOfCenters and h intersect
  Eigen::Vector2d hIntersect = h.cross(lineOfCenters).hnormalized();

  double xmax, xmin, ymax, ymin;
  xmax = (centerA(0)>centerB(0)) ? centerA(0) : centerB(0);
  xmin = (centerA(0)<centerB(0)) ? centerA(0) : centerB(0);
  ymax = (centerA(1)>centerB(1)) ? centerA(1) : centerB(1);
  ymin = (centerA(1)<centerB(1)) ? centerA(1) : centerB(1);
  bool g_fits, h_fits;
  
  g_fits = gIntersect(0)>xmin && gIntersect(0)<xmax && gIntersect(1)>ymin && gIntersect(1)<ymax;
  h_fits = hIntersect(0)>xmin && hIntersect(0)<xmax && hIntersect(1)>ymin && hIntersect(1)<ymax;
  if (!(g_fits ^ h_fits)) {
    std::cerr << "G|H -> " << g_fits << "|" << h_fits << std::endl;
    std::cerr << "Both or neither g_fits and h_fits are valid." << std::endl;
    return false;
  }
  l = g_fits ? g : h;
  return true;
}

bool computeInvariant(const Eigen::Vector3d& line1, 
                      const Eigen::Vector3d& line2, 
                      const Eigen::Matrix3d& locus,
                      double& invariant) {
  // if(vectorContainsNaN(line1) || vectorContainsNaN(line2)) {
  //     std::cerr << "One of the lines contains a NaN." << std::endl;
  //     return false;
  // }
  Eigen::Matrix3d envelope = adjugate(locus);
  // the numerator can be negative
  Eigen::RowVector3d line1T = line1.transpose();
  Eigen::RowVector3d line2T = line2.transpose();
  Eigen::MatrixXd numerator = line1T * envelope * line2;
  double numeratord = numerator.value();
  // the denominator is quadratic in both terms, so never negative
  Eigen::MatrixXd l1tel1 = line1T*envelope*line1;
  Eigen::MatrixXd l2tel2 = line2T*envelope*line2;
  double l1tel1d = l1tel1.value();
  double l2tel2f = l2tel2.value();
  double denominator = sqrt(l1tel1d*l2tel2f);
  // calculate the invariant, eqns 95-97
  invariant = acosh(abs(numeratord)/denominator);
  return true;
}

bool computeCraterTriadInvariants(const Conic& A, const Conic& B, const Conic& C,
                                  std::array<double, NONCOPLANAR_INVARIANTS>& invariants) {
  std::tuple<Eigen::Vector3d, Eigen::Vector3d> gh_ij, gh_jk, gh_ki;
  Eigen::Vector3d lij, ljk, lki;
  double invA, invB, invC;
  Eigen::Matrix3d locusA, locusB, locusC;
  Eigen::Vector2d centerA, centerB, centerC;
  locusA = A.getLocus();
  locusB = B.getLocus();
  locusC = C.getLocus();
  centerA = A.getCenter();
  centerB = B.getCenter();
  centerC = C.getCenter();

  if(!A.chooseIntersection(B, lij)) {
    std::cerr<<"IntersectionLines error ij\n";
    return false;
  }
  if(!B.chooseIntersection(C, ljk)) {
    std::cerr<<"IntersectionLines error jk\n";
    return false;
  }
  if(!C.chooseIntersection(A, lki)) {
    std::cerr<<"IntersectionLines error ki\n";
    return false;
  }
  // Invariants
  if(!computeInvariant(lij, lki, locusA, invA)) {
    std::cerr<<"computeInvariant error A\n";
    return false;
  }
  if(!computeInvariant(ljk, lij, locusB, invB)) {
    std::cerr<<"computeInvariant error B\n";
    return false;
  }
  if(!computeInvariant(lki, ljk, locusC, invC)) {
    std::cerr<<"computeInvariant error C\n";
    return false;
  }
  invariants = {invA, invB, invC};
  // invariants.push_back(invA);
  // invariants.push_back(invB);
  // invariants.push_back(invC);
  return true;
}

Eigen::Matrix3d canonical(const Eigen::Matrix3d& image_conic) {
  double focal_length = 1;
  // Applying the concept of the canonical frame from Kanatani:
  // "3D interpretation of conics and orthogonality"
  double A = image_conic(0, 0);
  double B = image_conic(0, 1);
  double C = image_conic(1, 1);
  // D = image_conic[0, 2]/focal_length
  // E = image_conic[1, 2]/focal_length
  // F = image_conic[2, 2]/focal_length**2
  // if AC-B^2 > 0:
  //   - if A+C > 0, ellipse
  //   - if A+C < 0, imaginary conic
  // if AC-B^2  < 0, hyperbola
  // if AC-B^2 == 0, parabola
  double ACmB2 = A*C - std::pow(B,2);
  // this value must not be zero for an ellipse
  assert(std::abs(ACmB2)>1e-40);
  double ApC = A+C;
  assert( ApC > 0);
  double lambda1 = (ApC - std::sqrt(std::pow(ApC,2) - 4*ACmB2))/2;
  double lambda2 = (ApC + std::sqrt(std::pow(ApC,2) - 4*ACmB2))/2;
  double mu = std::pow(focal_length,2)/ACmB2;
  // to be an ellipse, mu*lambda_i must be > 0
  assert(mu*lambda1 > 0 and mu*lambda2 > 0);

  double a = std::sqrt(mu/lambda1);
  double b = std::sqrt(mu/lambda2);
  double a2 = std::pow(a,2);
  double b2 = std::pow(b,2);
  double A_prime = b2;
  double C_prime = a2;
  double F_prime = -a2*b2;
  Eigen::Matrix3d canonical = Eigen::Matrix3d::Identity();
  canonical.diagonal() = Eigen::Vector3d(A_prime, C_prime, F_prime);
  return canonical;
}

Conic canonical(const Conic& conic) {
  Eigen::Matrix3d mtx = canonical(conic.getLocus());
  Conic canonical_conic(mtx);
  return canonical_conic;
}

// def calc_t1(lambda_i, conic_locus):
//     b = conic_locus[1, 1]
//     g = conic_locus[0, 2]
//     f = conic_locus[1, 2]
//     h = conic_locus[0, 1]
//     return (b-lambda_i)*g - f*h


// def calc_t2(lambda_i, conic_locus):
//     a = conic_locus[0, 0]
//     g = conic_locus[0, 2]
//     f = conic_locus[1, 2]
//     h = conic_locus[0, 1]
//     return (a-lambda_i)*f - g*h


// def calc_t3(lambda_i, conic_locus, t1, t2):
//     a = conic_locus[0, 0]
//     g = conic_locus[0, 2]
//     h = conic_locus[0, 1]
//     return -(a-lambda_i)*(t1/t2)/g - h/g


// def calc_t_values(lambda_i, conic_locus):
//         t1 = calc_t1(lambda_i, conic_locus)
//         t2 = calc_t2(lambda_i, conic_locus)
//         t3 = calc_t3(lambda_i, conic_locus, t1, t2)
//         return t1, t2, t3


// def calc_camera_xform(conic_locus, lambda_arr, n, m, l):
//     l_arr = np.empty((3,))
//     m_arr = np.empty((3,))
//     n_arr = np.empty((3,))
//     for index, lambda_i in enumerate(lambda_arr):
//         // TODO: combine calculating t_i
//         t1, t2, t3 = calc_t_values(lambda_i, conic_locus)
//         m_arr[index] = 1/np.sqrt(1+(t1/t2)**2+t3**2)
//         l_arr[index] = t1/t2*m_arr[index]
//         n_arr[index] = t3*m_arr[index]
//         pass
//     xform = np.row_stack((l_arr, m_arr, n_arr))
//     return xform

// Eigen::Matrix3d get_canonical_transform(const Eigen::Matrix3d& image_conic_locus, const double radius) {
//   // Three-dimensional ... machine vision
//   // eqns 30, 31, 33
//   Conic conic(image_conic_locus);

//   double Rmax = conic.getSemiMajorAxis();
//   double Rmin = conic.getSemiMinorAxis();
//   double e=1;
//   double e2 = std::pow(e, 2);
//   double Rmax2 = std::pow(Rmax, 2);
//   double Rmin2 = std::pow(Rmin, 2);

//   double alpha_angle = std::asin(std::sqrt((e2/Rmin2 - e2/Rmax2)/(e2/Rmin2 + 1)));
//   std::cout << "Angle: " << rad2deg(alpha_angle) << std::endl;
//   double dist = radius*(e/Rmin*std::cos(alpha_angle));
//   std::cout << "Distance: " << dist << std::endl;
  
//   Eigen::Vector2d n;
//   Eigen::Vector2d m;
//   Eigen::Vector2d l;
//   Eigen::Vector3d eigval = conic.getLocus().eigenvalues();
//   // eigenvals = np.diag(eigval)
//   // print(image_conic)
//   // print(eigenvals)
//   double lambda1 = eigval(0);
//   double lambda2 = eigval(1);
//   double lambda3 = eigval(2);
//   if(lambda1 < lambda2) {
//     n(0) = std::sqrt((lambda1-lambda3)/(lambda2-lambda3));
//     n(1) = n(0);
//     m(0) = std::sqrt((lambda2-lambda1)/(lambda2-lambda3));
//     m(1) = -m(0);
//     l(0) = 0;
//     l(1) = l(0);
//   }
//   else if(lambda1 > lambda2) {
//     n[0] = std::sqrt((lambda2-lambda3)/(lambda1-lambda3));
//     n[1] = n[0];
//     m[0] = 0;
//     m[1] = m[0];
//     l[0] = std::sqrt((lambda1-lambda2)/(lambda1-lambda3));
//     l[1] = -l[0];
//   }
//   else if(lambda1 == lambda2) { // right circular cone
//     n(0) = 1;
//     n(1) = n(0);
//     m(0) = 0;
//     m(1) = m(0);
//     l(0) = 0;
//     l(1) = l(0);
//   }
//   Eigen::Vector3d lambda_arr(lambda1, lambda2, lambda3);
//   xform = calc_camera_xform(image_conic, lambda_arr, n, m, l)
//   // print(xform)
//   // resulting_conic = xform.T@image_conic@xform
//   print(xform.T@image_conic@xform)
//   print(xform@image_conic@xform.T)
//   return xform
// }

} // namespace
