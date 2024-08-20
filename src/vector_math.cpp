#include "vector_math.h"

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
      //  element which are not in given row and
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

Eigen::MatrixXd getCofactorMatrix(const Eigen::MatrixXd& matrix) {
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

Eigen::MatrixXd getAdjugateMatrix(const Eigen::MatrixXd& matrix) {
  // TODO: ensure that row and column counts are equal
  Eigen::MatrixXd mtx_adj = getCofactorMatrix(matrix);
  // Matrix adjugate is the transpose of the cofactor matrix
  return mtx_adj.transpose();
}

Eigen::Matrix3d get3x3SymmetricAdjugateMatrix(const Eigen::Matrix3d& mtx) {
  Eigen::Matrix3d mtx_adj = Eigen::Matrix3d(3, 3);

  // get elements of A
  const auto a_11 = mtx.coeff(0, 0);
  const auto a_12 = mtx.coeff(0, 1);
  const auto a_13 = mtx.coeff(0, 2);

  const auto a_22 = mtx.coeff(1, 1);
  const auto a_23 = mtx.coeff(1, 2);

  const auto a_33 = mtx.coeff(2, 2);

  // Compute entries of cofactor matrix which, in a symmetric matrix, 
  // are the entries of the adjugate
  // entries Aij == Aji
  const auto mtx_adj11 =  a_22*a_33 - a_23*a_23;
  const auto mtx_adj12 = -a_12*a_33 + a_23*a_13;
  const auto mtx_adj13 =  a_12*a_23 - a_22*a_13;

  const auto mtx_adj21 = mtx_adj12;
  const auto mtx_adj22 = a_11*a_33 - a_13*a_13;
  const auto mtx_adj23 = -a_11*a_23 + a_12*a_13;

  const auto mtx_adj31 = mtx_adj13;
  const auto mtx_adj32 = mtx_adj23;
  const auto mtx_adj33 = a_11*a_22 - a_12*a_12;
  mtx_adj << mtx_adj11, mtx_adj12, mtx_adj13,
              mtx_adj21, mtx_adj22, mtx_adj23, 
              mtx_adj31, mtx_adj32, mtx_adj33;
  return mtx_adj;
}

double vdot(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2) {
  return point1.dot(point2);
}

double angularPseudoDistance(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2) {
  return vdot(point1, point2);
}

double angularDistance(const Eigen::Vector3d& point1, const Eigen::Vector3d& point2) {
  return acos(angularPseudoDistance(point1, point2));
}

// template <typename T>
// bool normalizeDeterminant(Eigen::DenseBase<T>& mtx) {
bool normalizeDeterminant(Eigen::Matrix4d& mtx) {
  // using approach from Matrix Cookbook
  // https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf
  // get the determinant
  uint ncol = mtx.cols();
  // uint nrow = mtx.rows();
  // T mtx_normalized(nrow, ncol);
  //get location of maximum
  // Eigen::Index maxRow, maxCol;
  double maxVal = 1/mtx.maxCoeff();
  mtx *= maxVal;
  double det_mtx = mtx.determinant();
  if(abs(det_mtx)<1e-20) {
    // BOOST_LOG_TRIVIAL(warning) << "Matrix is singular.";
    std::cerr << "Matrix is singular/nearly singular." << std::endl;
    mtx *= maxVal;
    return false;
  }
  // we want the determinant of A to be 1
  // 1 = det(c*A) = d^n*det(A), where n=3 for 3x3 matrix
  // so solve the equation d^3 - 1/det(A) = 0
  double d = pow(abs(1./det_mtx), 1./ncol);
  // d*mtx should now have a determinant of +1
  mtx *= d;
  bool sign = signbit(mtx.determinant()); // true if negative
  // assert(abs(det(A_norm)-1)<sqrt(eps))
  mtx *= sign ? -1.0 : 1.0;
  return true;
}

Eigen::Matrix3d crossMatrix(const Eigen::Vector3d& v) {
  Eigen::Matrix3d v_cross(3, 3);
  v_cross << 0, -v(2), v(1),
             v(2), 0, -v(0),
            -v(1), v(0), 0;
  return v_cross;
}

void normalizeVector(Eigen::Vector3d& vec) {
  vec /= vec.norm();
}

void normalizeVector(const Eigen::Vector3d& inVec, Eigen::Vector3d& outVec) {
  outVec = inVec;
  normalizeVector(outVec);
}

Eigen::Vector3d normalizeVector(const Eigen::Vector3d& inVec) {
  return inVec / inVec.norm();
}

Eigen::Vector3d GetNorthPoleUnitVector() {
  return Eigen::Vector3d::UnitZ();
}

void GetNorthPoleUnitVector(Eigen::Vector3d& north_pole) {
  north_pole = GetNorthPoleUnitVector();
}
