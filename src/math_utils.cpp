
#include "math_utils.h"
#include <iostream>
#include <vector>

// template <typename T>
// void operator /(Eigen::Vector3d& vec, T divisor) {
//   size_t siz = vec.size();
//   // for(auto& item : vec) {
//   for(size_t i=0; i<siz; i++) {
//     vec(i) /= divisor;
//     // item /= divisor;
//   }
// }

// template <typename T>
// void operator *(Eigen::Vector3d& vec, T scalar) {
//   size_t siz = vec.size();
//   // for(auto& item : vec) {
//   //   item *= scalar;
//   // }
//   for(size_t i=0; i<siz; i++) {
//     vec(i) *= scalar;
//   }
// }

double getCofactor(const Eigen::MatrixXd& matrix, size_t cf_row, size_t cf_col) {
  size_t nrow = matrix.rows();
  size_t ncol = matrix.cols();
  size_t i = 0, j = 0;
  Eigen::MatrixXd cf_temp(nrow-1, ncol-1);

  // Looping for each element of the matrix
  for (size_t irow = 0; irow < nrow; irow++) {
    for (size_t icol = 0; icol < ncol; icol++) {
      //  Copying into temporary matrix only those
      //  elements which are not in given row and
      //  column
      if (irow != cf_row && icol != cf_col) {
        cf_temp(i, j++) = matrix(irow, icol);

        // Row is filled, so increase row index and reset col index
        if (j == nrow - 1) {
          j = 0;
          i++;
        }
      }
    }
  }
  return cf_temp.determinant();
}

Eigen::MatrixXd cofactor(const Eigen::MatrixXd& matrix) {
  size_t nrow = matrix.rows();
  size_t ncol = matrix.cols();
  Eigen::MatrixXd cofactor_matrix(nrow, ncol);
  int cofactor_sign;
  for(size_t row=0;row<nrow; row++) {
    for(size_t col=0; col<ncol; col++) {
      cofactor_sign = (row+col)%2==0 ? 1 : -1;
      cofactor_matrix(row, col) = cofactor_sign*getCofactor(matrix, row, col);
    }
  }
  return cofactor_matrix;
}

Eigen::MatrixXd adjugate(const Eigen::MatrixXd& matrix) {
  // TODO: ensure that row and column counts are equal
  // std::cout << __func__ << " Matrix:\n" << matrix << std::endl;
  if(matrix.rows() != matrix.cols()) {
    throw std::runtime_error("Matrix is not square");
  }
  if(matrix.rows() == 3 && matrix.cols() == 3 && matrix.isApprox(matrix.transpose())) {
    return symmetricAdjugate(matrix);
  }
  Eigen::MatrixXd mtx_adj = cofactor(matrix);

  // std::cout << __func__ << " Adjugate:\n" << mtx_adj << std::endl;
  if(mtx_adj.hasNaN()) {
    throw std::runtime_error("Matrix contains NaN values.");
  }
  // Just remember: 
  //  "The adjugate of my adjugate is a factor of the determinant to a power."
  // Matrix adjugate is the transpose of the cofactor matrix
  // Just as a fun fact, the adj(adj(mtx)) = det(mtx)^(n-2) 
  // where n is the number of dimensions (e.g., 4x4 => 4 dim)
  // Unfortunately, some of the work with quadrics uses disk
  // quadrics, which are rank deficient (i.e., det=0), 
  // and this fun fact does not hold
  // Another check for adjugates is ensure that the matrix times its
  // adjugate is the identity matrix scaled by the det of the matrix
  mtx_adj = (mtx_adj.cwiseAbs().array() < 1e-100).select(0, mtx_adj);

  return mtx_adj.transpose();
}

Eigen::Matrix3d symmetricAdjugate(const Eigen::Matrix3d& mtx) {
  if(mtx.hasNaN()) {
    throw std::runtime_error("Matrix contains NaN values.");
  }
  Eigen::Matrix3d mtx_adj = Eigen::Matrix3d(3, 3);

  // get elements of A
  const double a_11 = mtx.coeff(0, 0);
  const double a_12 = mtx.coeff(0, 1);
  const double a_13 = mtx.coeff(0, 2);

  const double a_22 = mtx.coeff(1, 1);
  const double a_23 = mtx.coeff(1, 2);

  const double a_33 = mtx.coeff(2, 2);

  // Compute entries of cofactor matrix which, in a symmetric matrix, 
  // are the entries of the adjugate
  // entries Aij == Aji
  const double mtx_adj11 =  a_22*a_33 - a_23*a_23;
  const double mtx_adj12 = -a_12*a_33 + a_23*a_13;
  const double mtx_adj13 =  a_12*a_23 - a_22*a_13;

  const double mtx_adj21 = mtx_adj12;
  const double mtx_adj22 = a_11*a_33 - a_13*a_13;
  const double mtx_adj23 = -a_11*a_23 + a_12*a_13;

  const double mtx_adj31 = mtx_adj13;
  const double mtx_adj32 = mtx_adj23;
  const double mtx_adj33 = a_11*a_22 - a_12*a_12;
  mtx_adj <<  mtx_adj11, mtx_adj12, mtx_adj13,
              mtx_adj21, mtx_adj22, mtx_adj23, 
              mtx_adj31, mtx_adj32, mtx_adj33;
  return mtx_adj;
}

double vdot(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2) {
  return point1.dot(point2);
}

double getPseudoAngleBetweenVectors(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2) {
  return vdot(point1, point2);
}

double getAngleBetweenVectors(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2) {
  return acos(getPseudoAngleBetweenVectors(point1, point2));
}


Eigen::Vector3d getAxisNormalToVectors(const Eigen::Vector3d& vec1, const Eigen::Vector3d& vec2) {
  if(vec1.isApprox(vec2)) {
    throw std::runtime_error("Vectors are nearly parallel; invalid axis");
  }
  return vec1.cross(vec2).normalized();
}

Eigen::Matrix3d crossMatrix(const Eigen::Vector3d& v) {
  Eigen::Matrix3d v_cross(3, 3);
  v_cross << 0, -v(2), v(1),
             v(2), 0, -v(0),
            -v(1), v(0), 0;
  return v_cross;
}

Eigen::Vector3d getNorthPoleUnitVector() {
  return Eigen::Vector3d::UnitZ();
}

void GetNorthPoleUnitVector(Eigen::Vector3d& north_pole) {
  north_pole = getNorthPoleUnitVector();
}

Eigen::Matrix3d getENUFrame(const Eigen::Vector3d& surface_point) {
  Eigen::Vector3d u_north_pole(3), u_surface_point(3), ui(3), ei(3), ni(3);
  Eigen::Matrix3d enu_frame(3, 3);
  // double eps = 1e-6;
  u_north_pole = getNorthPoleUnitVector();
  u_surface_point = surface_point.normalized();
  // TODO: should latitudes near the pole be errors?
  // TODO: we could also replace acos() call with dot product minus 1 is approx zero
  if(acos(u_north_pole.dot(u_surface_point))<EPS) {
    // std::cerr << "Surface point is near the North Pole.\n";
    ei << 1, 0, 0;
    ni << 0, 1, 0;
    ui = u_north_pole;
    throw std::runtime_error("Surface point is near the North Pole.");
  }  
  // TODO: we could also replace acos() call with dot product minus 1 is approx zero OR
  // alternatively, u_north_pole.dot(u_surface_point) + 1 is approx zero
  else if(acos(-u_north_pole.dot(u_surface_point))<EPS) {
    // std::cerr << "Surface point is near the South Pole.\n";
    ei << -1, 0, 0;
    ni << 0, -1, 0;
    ui = -u_north_pole;
    throw std::runtime_error("Surface point is near the South Pole.");
  }
  else {
    ui = u_surface_point;
    ei = u_north_pole.cross(ui);
    ni = ui.cross(ei);
  }
  enu_frame << ei.normalized(), ni.normalized(), ui.normalized();
  return enu_frame;
}
Eigen::Matrix3d getENUFrame(const double lat, const double lon) {
  if(std::abs(lat-90) <EPS) {
    throw std::runtime_error("Surface point is near the North Pole.");
  }
  else if(std::abs(lat+90) <EPS) {
    throw std::runtime_error("Surface point is near the South Pole.");
  }
  const double lambda = deg2rad(lon);
  const double phi    = deg2rad(lat);
  Eigen::Matrix3d mtx;
  const double slambda = std::sin(lambda);
  const double sphi    = std::sin(phi);
  const double clambda = std::cos(lambda);
  const double cphi    = std::cos(phi);
  mtx << -slambda, -clambda*sphi, clambda*cphi,
          clambda, -slambda*sphi, slambda*cphi,
                0,          cphi,         sphi;
  return mtx;                
}

// Eigen::Matrix3d getENUFrame(const double lat, const double lon) {
//   Eigen::Vector3d u_north_pole(3), ui(3), ei(3), ni(3);
//   Eigen::Matrix3d enu_frame(3, 3);
//   u_north_pole = getNorthPoleUnitVector();
//   // TODO: should latitudes near the pole be errors?
//   if(std::abs(lat-90) <EPS) {
//     // std::cerr << "Surface point is near the North Pole.\n";
//     ei << 1, 0, 0;
//     ni << 0, 1, 0;
//     ui = u_north_pole;
//     throw std::runtime_error("Surface point is near the North Pole.");
//   }
//   else if(std::abs(lat+90) <EPS) {
//     // std::cerr << "Surface point is near the South Pole.\n";
//     ei << -1, 0, 0;
//     ni << 0, -1, 0;
//     ui = -u_north_pole;
//     throw std::runtime_error("Surface point is near the South Pole.");
//   }
//   else {
//     ui = latlon2bearing(lat, lon);
//     ei = u_north_pole.cross(ui);
//     ni = ui.cross(ei);
//   }
//   enu_frame << ei.normalized(), ni.normalized(), ui.normalized();
//   return enu_frame;
// }

Eigen::Quaterniond eulerToQuaternion(const double roll,
                                     const double pitch,
                                     const double yaw) {
  Eigen::AngleAxisd rollAngle(roll, Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd pitchAngle(pitch, Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd yawAngle(yaw, Eigen::Vector3d::UnitZ());

  Eigen::Quaterniond q = yawAngle * pitchAngle * rollAngle;
  return q;
}

void eulerToDCM(const double roll,
                const double pitch,
                const double yaw,
                Eigen::Matrix3d& dcm) {
  Eigen::AngleAxisd rollAngle(roll, Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd pitchAngle(pitch, Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd yawAngle(yaw, Eigen::Vector3d::UnitZ());

  dcm = yawAngle * pitchAngle * rollAngle;
}

// Templated version is in header
void convertEigenVectorToVector(const Eigen::Vector3d& eig, std::vector<double>& vec) {
  try {
    Eigen::Vector3d::Map(&vec[0], eig.size()) = eig;
  }
  catch (const std::exception& e) {
    std::cerr << "Could not convert " << eig.transpose() << " to ";
    for (const auto& elem : vec) {
      std::cerr << elem << ", ";
    }
    std::cout << std::endl;
    throw std::runtime_error("Cannot convert eigen3 vector to std::vector.");
  }
}

// Templated version is in header
bool vectorContainsNaN(const Eigen::Vector3d& vec) {
  return vec.hasNaN();
}
