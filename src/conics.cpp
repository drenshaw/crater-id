#include <iostream>
// #include <stdexcept>
#include <tuple>
#include <array>
#include <math.h>
#include <eigen3/Eigen/Dense>

#include "conics.h"
#include "math_utils.h"

int Conic::next_id = 0;

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
Conic::Conic(const std::vector<double>& vec) : 
  semimajor_axis_{vec.at(0)},
  semiminor_axis_{vec.at(1)},
  x_center_      {vec.at(2)},
  y_center_      {vec.at(3)},
  angle_         {vec.at(4)} {
    setID();
}

void Conic::setID() {
  id_ = next_id++;
}

int Conic::GetID() const {
  return id_;
}

bool Conic::operator==(const Conic& other_conic) const {
  return  almost_equal(this->semimajor_axis_, other_conic.GetSemiMajorAxis()) &&
          almost_equal(this->semiminor_axis_, other_conic.GetSemiMinorAxis()) &&
          almost_equal(this->x_center_, other_conic.GetCenterX()) &&
          almost_equal(this->y_center_, other_conic.GetCenterY()) &&
          almost_equal(this->angle_, other_conic.GetAngle());
}

bool Conic::operator!=(const Conic& other_conic) const {
  return  !almost_equal(semimajor_axis_, other_conic.semimajor_axis_) ||
          !almost_equal(semiminor_axis_, other_conic.semiminor_axis_) ||
          !almost_equal(x_center_, other_conic.x_center_) ||
          !almost_equal(y_center_, other_conic.y_center_) ||
          !almost_equal(angle_, other_conic.angle_);
}

std::ostream& operator<<(std::ostream& os, const Conic& conic) {
    return os 
        << "Conic-> \tID: '" << conic.id_ << "'"
        << "\n\tSemi axes: " 
        << conic.semimajor_axis_ << " | "
        << conic.semiminor_axis_
        << "\n\tCenter: (" << conic.x_center_ << " , "
        << conic.y_center_ << ") "
        << "\n\tAngle (deg): " << rad2deg(conic.angle_) << " ";
}

void Conic::SetGeometricParameters(const std::array<double, GEOMETRIC_PARAM>& geom_vec) {
  SetGeometricParameters( geom_vec.at(0), 
                          geom_vec.at(1), 
                          geom_vec.at(2), 
                          geom_vec.at(3), 
                          geom_vec.at(4));
}

void Conic::SetGeometricParameters(const std::vector<double>& geom_vec) {
  SetGeometricParameters( geom_vec.at(0), 
                          geom_vec.at(1), 
                          geom_vec.at(2), 
                          geom_vec.at(3), 
                          geom_vec.at(4));
}

void Conic::SetGeometricParameters( const double semimajor_axis, 
                                    const double semiminor_axis, 
                                    const double x_center, 
                                    const double y_center, 
                                    const double angle) {
  semimajor_axis_ = semimajor_axis;
  semiminor_axis_ = semiminor_axis;
  x_center_ = x_center;
  y_center_ = y_center;
  angle_ = angle;
}

void Conic::SetSemimajorAxis(const double semimajor_axis) {
  semimajor_axis_ = semimajor_axis;
}

void Conic::SetSemiminorAxis(const double semiminor_axis) {
  semiminor_axis_ = semiminor_axis;
}

void Conic::SetCenterX(const double x_center) {
  x_center_ = x_center;
}

void Conic::SetCenterY(const double y_center) {
  y_center_ = y_center;
}

void Conic::SetAngle(const double angle) {
  angle_ = angle;
}


void Conic::SetImplicitParameters(const std::array<double, IMPLICIT_PARAM>& impl_params) {
  std::array<double, GEOMETRIC_PARAM> geom_params = Implicit2Geom(impl_params);
  SetGeometricParameters(geom_params);
}

void Conic::SetLocus(const Eigen::Matrix3d& locus) {
  std::array<double, GEOMETRIC_PARAM> geom_params = Locus2Geom(locus);
  SetGeometricParameters(geom_params);
}

void Conic::NormalizeImplicitParameters(std::array<double, IMPLICIT_PARAM>& impl_params) const {
  // double vecNormRecip = 1/vectorNorm(impl_params);
  double vecNormRecip = 1/impl_params.at(5);
  std::transform(impl_params.begin(), impl_params.end(), impl_params.begin(),
               std::bind(std::multiplies<double>(), std::placeholders::_1, vecNormRecip));
}

// TODO: possibly create template or call to same function for std::array above
void Conic::NormalizeImplicitParameters(std::vector<double>& impl_params) const {
  // double vecNormRecip = 1/vectorNorm(impl_params);
  double vecNormRecip = 1/impl_params.at(5);
  std::transform(impl_params.begin(), impl_params.end(), impl_params.begin(),
               std::bind(std::multiplies<double>(), std::placeholders::_1, vecNormRecip));
}

std::array<double, GEOMETRIC_PARAM> Conic::GetGeom() const {
  return {semimajor_axis_, semiminor_axis_, x_center_, y_center_, angle_};
}

std::array<double, IMPLICIT_PARAM> Conic::GetImplicit() const {
  return Geom2Implicit();
}

Eigen::Matrix3d Conic::GetLocus() const {
  return Geom2Locus();
}

Eigen::Matrix3d Conic::GetEnvelope() const {
  return getAdjugateMatrix(Geom2Locus());
}

std::array<double, IMPLICIT_PARAM> Conic::Locus2Implicit(const Eigen::Matrix3d& locus) const {
  return {
    locus.coeff(0,0),
    2*locus.coeff(0,1),
    locus.coeff(1,1),
    2*locus.coeff(0,2),
    2*locus.coeff(1,2),
    locus.coeff(2,2),
  };
}

std::array<double, GEOMETRIC_PARAM> Conic::Implicit2Geom(const std::array<double, IMPLICIT_PARAM>& impl_params) const {
  double numerator, denominator_a, denominator_b;
  double B2_minus_4AC, xc, yc, phi;
  double semimajor_axis, semiminor_axis, amc2, b2;
  double A = impl_params.at(0);
  double B = impl_params.at(1);
  double C = impl_params.at(2);
  double D = impl_params.at(3);
  double E = impl_params.at(4);
  double F = impl_params.at(5);
  // Compute discriminant for quadratic equation
  B2_minus_4AC = pow(B, 2) - 4*A*C;
  // Compute ellipse center (See Eq 4.16 in [Christian, 2010])
  xc = (2 * C * D - B * E) / B2_minus_4AC;
  yc = (2 * A * E - B * D) / B2_minus_4AC;
  // Compute ellipse semimajor axis (a) and seminor axis (b)
  // (See Eqs 4.17 and 4.18 in [Christian, 2010])
  numerator = 2*(A*E*E + C*D*D - B*D*E + F*B2_minus_4AC);
  amc2 = pow(A - C, 2);
  b2 = pow(B, 2);
  denominator_a = B2_minus_4AC*( sqrt(amc2 + b2) - A - C);
  denominator_b = B2_minus_4AC*(-sqrt(amc2 + b2) - A - C);
  semimajor_axis = sqrt(numerator / denominator_a);
  semiminor_axis = sqrt(numerator / denominator_b);
  // Compute angle from the x axis to semimajor axis direction
  // (See Eq 4.19 in [Christian, 2010])
  if(B == 0) {
    phi = 0.0;
  } else {
    phi = 0.5 * atan2(B, A - C);
  }
  if(A>C) {
    phi += M_PI / 2;
  }
  return {semimajor_axis, semiminor_axis, xc, yc, phi};
}

Eigen::Matrix3d Conic::Geom2Locus() const {
  const std::array<double, IMPLICIT_PARAM> impl_params = Geom2Implicit();
  // return normalize_determinant(Implicit2Locus(impl_params))
  return Implicit2Locus(impl_params);
 }

Eigen::Matrix3d Conic::Implicit2Locus(const std::array<double, IMPLICIT_PARAM>& impl_params) const {
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

Eigen::Matrix3d Conic::Implicit2Locus() const {
  return Implicit2Locus(Geom2Implicit());
}

 std::array<double, GEOMETRIC_PARAM> Conic::Locus2Geom(const Eigen::Matrix3d& locus) const {
  const std::array<double, IMPLICIT_PARAM> impl_params = Locus2Implicit(locus);
  return Implicit2Geom(impl_params);
}


std::array<double, IMPLICIT_PARAM> Conic::Geom2Implicit() const {
  // unpack geometric parameters
  double a   = semimajor_axis_;
  double b   = semiminor_axis_;
  double xc  = x_center_;
  double yc  = y_center_;
  double phi = angle_;

  // perform some computations beforehand
  double a2 = pow(a, 2);
  double b2 = pow(b, 2);
  double sin_phi = sin(deg2rad(phi));
  double cos_phi = cos(deg2rad(phi));
  double xc_2 = pow(xc, 2);
  double yc_2 = pow(yc, 2);
  double sin_phi_2 = pow(sin_phi, 2);
  double cos_phi_2 = pow(cos_phi, 2);

  std::array<double, IMPLICIT_PARAM> coeff;

  coeff.at(0) =  a2*sin_phi_2 + b2*cos_phi_2;
  coeff.at(1) =  2*(b2-a2)*sin_phi*cos_phi;
  coeff.at(2) =  a2*cos_phi_2 + b2*sin_phi_2;
  coeff.at(3) = -2*coeff[0]*xc - coeff[1]*yc;
  coeff.at(4) = -coeff[1]*xc - 2*coeff[2]*yc;
  coeff.at(5) =  coeff[0]*xc_2 + coeff[1]*xc*yc + coeff[2]*yc_2 - a2*b2;

  NormalizeImplicitParameters(coeff);
  return coeff;
}

bool Conic::conicIntersectionLines(const Eigen::Matrix3d& Aj, 
                                   std::tuple<Eigen::Vector3d, Eigen::Vector3d>& gh) const {
  Eigen::Matrix3d Ai = GetLocus();
  return invariants::intersectionLines(Ai, Aj, gh);
}

bool Conic::conicIntersectionLines(const Conic& conicB, 
                                   std::tuple<Eigen::Vector3d, Eigen::Vector3d>& gh) const {
  Eigen::Matrix3d Ai = GetLocus();
  Eigen::Matrix3d Aj = conicB.GetLocus();
  return invariants::intersectionLines(Ai, Aj, gh);
}

double Conic::GetCenterX() const {
  return x_center_;
}

double Conic::GetCenterY() const {
  return y_center_;
}

Eigen::Vector2d Conic::GetCenter() const {
  Eigen::Vector2d center;
  center << x_center_, y_center_;
  return center;
}

void Conic::GetCenter(cv::Point& center) const {
  center.x = x_center_;
  center.y = y_center_;
}

void Conic::GetCenter(Eigen::Vector2d& center) const {
  center << x_center_, y_center_;
}

double Conic::GetSemiMajorAxis() const {
  return semimajor_axis_;
}

double Conic::GetSemiMinorAxis() const {
  return semiminor_axis_;
}

cv::Size Conic::GetSemiAxes() const {
  cv::Size axes(semimajor_axis_, semiminor_axis_);
  return axes;
}

void Conic::GetSemiAxes(Eigen::Vector2d& semiaxes) const {
  semiaxes << semimajor_axis_, semiminor_axis_;
}

void Conic::GetSemiAxes(cv::Point& semiaxes) const {
  semiaxes.x = semimajor_axis_;
  semiaxes.y = semiminor_axis_;
}

void Conic::GetSemiAxes(cv::Size& semiaxes) const {
  semiaxes.width  = semimajor_axis_;
  semiaxes.height = semiminor_axis_;
}

// GetSize is an alias for GetSemiAxes
cv::Size Conic::GetSize() const {
  return GetSemiAxes();
}

void Conic::GetSize(cv::Size& semiaxes) const {
  return GetSemiAxes(semiaxes);
}

double Conic::GetAngle() const {
  return angle_;
}

bool Conic::chooseConicIntersection(const Conic& other, Eigen::Vector3d& l) const {
  std::tuple<Eigen::Vector3d, Eigen::Vector3d> gh;
  if(!invariants::intersectionLines(GetLocus(), other.GetLocus(), gh)) {
    std::cerr<<"IntersectionLines error ij\n";
    return false;
  }
  if(!invariants::chooseIntersection(gh, GetCenter(), other.GetCenter(), l)) {
    std::cerr<<"ChooseIntersection error jk\n";
    return false;
  }

  return true;
}

// void convertEigenVectorToVector(const Eigen::Vector3d& eig, std::array<double, CONIC_DIM>& arr) {
//   Eigen::Vector3d::Map(&arr[0], eig.size()) = eig;
// }

// void convertEigenVectorToVector(const Eigen::Vector3d& eig, std::vector<double>& vec) {
//   Eigen::Vector3d::Map(&vec[0], eig.size()) = eig;
// }

// // Templated version of vectorContainsNaN in header
// bool vectorContainsNaN(const Eigen::Vector3d& eV) {
//   std::vector<double> vec;
//   convertEigenVectorToVector(eV, vec);
//   return vectorContainsNaN(vec);
// }

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
  Bij_star = getAdjugateMatrix(Bij);

  bkk_eig = Bij_star.diagonal();
  std::vector<double> bkk;
  convertEigenVectorToVector(bkk_eig, bkk);
  // eq 86: all diagonal entries of Bij_star are negative
  if ( std::any_of(bkk.begin(), bkk.end(), [&eps](double i){return i>eps;}) ) {
      // std::cerr << "bkk contains positive numbers: \n" << bkk_eig << std::endl;
      return false;
  }
  int min_idx = arg_min(bkk);

  // eq 87
  Eigen::Vector3d z = -Bij_star.row(min_idx)/sqrt(-bkk.at(min_idx));
  // Eigen::Vector3d z = -Bij_star.row(min_idx)/sqrt(-min_elem);
  Eigen::Matrix3d D = Bij + crossMatrix(z);
  Eigen::Matrix3d Dabs = D.cwiseAbs();
  
  Eigen::Matrix3d::Index maxRow,maxCol;
  Dabs.maxCoeff(&maxRow, &maxCol);

  // pg 44 of Christian, Derksen, Watkins [2020]
  // g is non-zero column of D, h is non-zero row of D
  Eigen::Vector3d g = D.col(maxCol);
  Eigen::Vector3d h = D.row(maxRow);
  
  if ( vectorContainsNaN(g) ) {
      std::cout << "g vector contains NaN's: " << g << std::endl;
      return false;
  }
  if ( vectorContainsNaN(h) ) {
      std::cout << "h vector contains NaN's: " << h << std::endl;
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
  Eigen::Matrix3d envelope = getAdjugateMatrix(locus);
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

// TODO: Can't use "const" on the conics here; figure out how to do that
bool computeCraterTriadInvariants(const Conic& A, const Conic& B, const Conic& C,
                                  std::array<double, NONCOPLANAR_INVARIANTS>& invariants) {
  std::tuple<Eigen::Vector3d, Eigen::Vector3d> gh_ij, gh_jk, gh_ki;
  Eigen::Vector3d lij, ljk, lki;
  double invA, invB, invC;
  Eigen::Matrix3d locusA, locusB, locusC;
  Eigen::Vector2d centerA, centerB, centerC;
  locusA = A.GetLocus();
  locusB = B.GetLocus();
  locusC = C.GetLocus();
  centerA = A.GetCenter();
  centerB = B.GetCenter();
  centerC = C.GetCenter();

  if(!A.chooseConicIntersection(B, lij)) {
    std::cerr<<"IntersectionLines error ij\n";
    return false;
  }
  if(!B.chooseConicIntersection(C, ljk)) {
    std::cerr<<"IntersectionLines error jk\n";
    return false;
  }
  if(!C.chooseConicIntersection(A, lki)) {
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
} // namespace
