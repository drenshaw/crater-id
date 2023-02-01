
#include "conics.h"

Conic::Conic(const double semimajor_axis, 
             const double semiminor_axis, 
             const double x_center, 
             const double y_center, 
             const double angle) {
    semimajor_axis_ = semimajor_axis;
    semiminor_axis_ = semiminor_axis;
    x_center_ = x_center;
    y_center_ = y_center;
    angle_ = angle;
    locus_ = Geom2Locus();
}
Conic::Conic(const std::tuple<double, double, double, double, double>& geom_tup) {
    semimajor_axis_ = std::get<0>(geom_tup);
    semiminor_axis_ = std::get<1>(geom_tup);
    x_center_       = std::get<2>(geom_tup);
    y_center_       = std::get<3>(geom_tup);
    angle_          = std::get<4>(geom_tup);
    locus_          = Geom2Locus();
}
Conic::Conic(const std::vector<double>& geom_vec) {
    setGeometricParameters(geom_vec);
}

void Conic::setGeometricParameters(const std::vector<double>& geom_vec) {
    semimajor_axis_ = geom_vec[0];
    semiminor_axis_ = geom_vec[1];
    x_center_       = geom_vec[2];
    y_center_       = geom_vec[3];
    angle_          = geom_vec[4];
    locus_          = Geom2Locus();
}

void Conic::setImplicitParameters(const std::vector<double>& impl_params) {
    std::vector<double> geom_params = Implicit2Geom(impl_params);
    setGeometricParameters(geom_params);
}

void Conic::setLocus(const Eigen::Matrix3d& locus) {
    std::vector<double> geom_params = Locus2Geom(locus);
    setGeometricParameters(geom_params);
}

void Conic::NormalizeImplicitParameters(std::vector<double>& impl_params) {
    auto vecNorm = vectorNorm(impl_params);
    // std::cout << "Vector norm: " << vecNorm << std::endl;
    for(auto& element : impl_params) {
        element /= vecNorm;
    }
}

std::vector<double> Conic::getGeom() {
    return {semimajor_axis_, semiminor_axis_, x_center_, y_center_, angle_};
}

std::vector<double> Conic::getImplicit() {
    return Geom2Implicit();
}

Eigen::Matrix3d Conic::getLocus() {
    return locus_;
}

Eigen::Matrix3d Conic::getEnvelope() {
    return getMatrixAdjugate(Geom2Locus());
}

std::vector<double> Conic::Locus2Implicit(const Eigen::Matrix3d& locus) {
    std::vector<double> coeff;
    coeff.reserve(6);
    coeff.push_back(locus.coeff(0,0));
    coeff.push_back(2*locus.coeff(0,1));
    coeff.push_back(locus.coeff(1,1));
    coeff.push_back(2*locus.coeff(0,2));
    coeff.push_back(2*locus.coeff(1,2));
    coeff.push_back(locus.coeff(2,2));
    return coeff;
}

std::vector<double> Conic::Implicit2Geom(const std::vector<double>& impl_params){
    double numerator, denominator_a, denominator_b;
    double B2_minus_4AC, xc, yc, phi;
    double semimajor_axis, semiminor_axis, amc2, b2;
    double A, B, C, D, E, F;
    A = impl_params[0];
    B = impl_params[1];
    C = impl_params[2];
    D = impl_params[3];
    E = impl_params[4];
    F = impl_params[5];
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
        phi = A>C?M_PI/2. : 0.0;
    } else {
        phi = (1 / 2) * atan2(B, A - C);
        if(A>C) {
            phi += M_PI / 2;
        }
    }
    return {semimajor_axis, semiminor_axis, xc, yc, phi};
}

Eigen::Matrix3d Conic::Geom2Locus() {
    const std::vector<double> impl_params = Geom2Implicit();
    // return normalize_determinant(Implicit2Locus(impl_params))
    return Implicit2Locus(impl_params);
 }

Eigen::Matrix3d Conic::Implicit2Locus(const std::vector<double>& impl_params) {
    Eigen::Matrix3d locus_mtx(3,3);
    double A, B, C, D, E, F;
    A = impl_params[0];
    B = impl_params[1];
    C = impl_params[2];
    D = impl_params[3];
    E = impl_params[4];
    F = impl_params[5];
    locus_mtx << A,  B/2, D/2,
                B/2,  C,  E/2,
                D/2, E/2,  F;

    return locus_mtx;
}

 std::vector<double> Conic::Locus2Geom(const Eigen::Matrix3d& locus) {
    const std::vector<double> impl_params = Locus2Implicit(locus);
    return Implicit2Geom(impl_params);
 }


std::vector<double> Conic::Geom2Implicit() {
    // unpack geometric parameters
    double a   = semimajor_axis_;
    double b   = semiminor_axis_;
    double xc  = x_center_;
    double yc  = y_center_;
    double phi = angle_;

    // perform some computations beforehand
    double a2 = pow(a, 2);
    double b2 = pow(b, 2);
    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    double xc_2 = pow(xc, 2);
    double yc_2 = pow(yc, 2);
    double sin_phi_2 = pow(sin_phi, 2);
    double cos_phi_2 = pow(cos_phi, 2);

    // Populuate each coefficient
    std::vector<double> coeff;
    coeff.reserve(6);
    coeff.push_back( a2*sin_phi_2 + b2*cos_phi_2);
    coeff.push_back( 2*(b2-a2)*sin_phi*cos_phi);
    coeff.push_back( a2*cos_phi_2 + b2*sin_phi_2);
    coeff.push_back(-2*coeff[0]*xc - coeff[1]*yc);
    coeff.push_back(-coeff[1]*xc - 2*coeff[2]*yc);
    coeff.push_back( coeff[0]*xc_2 + coeff[1]*xc*yc + coeff[2]*yc_2 - a2*b2);

//     # normalize coefficients
    NormalizeImplicitParameters(coeff);
    return coeff;
}

bool Conic::ConicIntersectionLines(const Eigen::Matrix3d& Aj, 
                                   std::tuple<Eigen::Vector3d, Eigen::Vector3d>& gh) {
    Eigen::Matrix3d Ai = getLocus();
    return IntersectionLines(Ai, Aj, gh);
}

bool Conic::ConicIntersectionLines(Conic& conicB, 
                                   std::tuple<Eigen::Vector3d, Eigen::Vector3d>& gh) {
    Eigen::Matrix3d Ai = getLocus();
    Eigen::Matrix3d Aj = getLocus();
    return IntersectionLines(Ai, Aj, gh);
}

Eigen::Vector2d Conic::getCenter() {
    Eigen::Vector2d center;
    center << x_center_, y_center_;
    return center;
}

Eigen::Vector2d Conic::getSemiAxes() {
    Eigen::Vector2d axes;
    axes << semimajor_axis_, semiminor_axis_;
    return axes;
}

Eigen::Vector3d getNorthPoleUnitVector() {
    return Eigen::Vector3d::UnitZ();
}

void getNorthPoleUnitVector(Eigen::Vector3d& north_pole) {
    north_pole = getNorthPoleUnitVector();
}

std::vector<double> convertEigenVectorToVector(const Eigen::Vector3d& eig) {
    std::vector<double> regVec;
    regVec.resize(eig.size());
    Eigen::Vector3d::Map(&regVec[0], eig.size()) = eig;
    return regVec;
}

template <typename T>
bool vectorContainsNaN(const std::vector<T> vec) {
    return std::any_of(vec.begin(), vec.end(), [](double i){return std::isnan(i);});
}

bool vectorContainsNaN(const Eigen::Vector3d& eV) {
    std::vector<double> vec = convertEigenVectorToVector(eV);
    return vectorContainsNaN(vec);
}

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

Eigen::MatrixXd getMatrixAdjugate(const Eigen::MatrixXd& matrix) {
    Eigen::MatrixXd mtx_adj = getCofactorMatrix(matrix);
    // Matrix adjugate is the transpose of the cofactor matrix
    return mtx_adj.transpose();
}

Eigen::Matrix3d get3x3SymmetricMatrixAdjugate(const Eigen::Matrix3d& mtx) {
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

bool IntersectConics(const Eigen::Matrix3d& Ai, 
                     const Eigen::Matrix3d& Aj, 
                     const double eig,
                     std::tuple<Eigen::Vector3d, Eigen::Vector3d>& gh) {
    double eps = 1e-16;
    Eigen::Matrix3d Bij(3, 3), Bij_star(3, 3);
    Eigen::Vector3d bkk_eig(3);
    Bij = eig*Ai + Aj;
    Bij_star = getMatrixAdjugate(Bij);

    bkk_eig = Bij_star.diagonal();
    std::vector<double> bkk = convertEigenVectorToVector(bkk_eig);
    // eq 86: all diagonal entries of Bij_star are negative
    if ( std::any_of(bkk.begin(), bkk.end(), [&eps](double i){return i>eps;}) ) {
        // std::cout << "bkk contains positive numbers: \n" << bkk_eig << std::endl;
        return false;
    }
    int min_idx = arg_min(bkk);

    // eq 87
    Eigen::Vector3d z = -Bij_star.row(min_idx)/sqrt(-bkk[min_idx]);
    Eigen::Matrix3d D = Bij + crossMatrix(z);
    Eigen::Matrix3d Dabs = D.cwiseAbs();
    
    Eigen::Matrix3d::Index maxRow,maxCol;
    Dabs.maxCoeff(&maxRow, &maxCol);

    // pg 44 of Christian, Derksen, Watkins [2020]
    // g is non-zero column of D, h is non-zero row of D
    Eigen::Vector3d g = D.col(maxCol);
    Eigen::Vector3d h = D.row(maxRow);
    
    std::vector<double> gVec = convertEigenVectorToVector(g);
    std::vector<double> hVec = convertEigenVectorToVector(h);
    if ( vectorContainsNaN(gVec) ) {
        std::cout << "g vector contains NaN's: " << g << std::endl;
        return false;
    }
    if ( vectorContainsNaN(hVec) ) {
        std::cout << "h vector contains NaN's: " << h << std::endl;
        return false;
    }
    gh = {g, h};
    return true;
}

bool IntersectionLines(const Eigen::Matrix3d& Ai, 
                            const Eigen::Matrix3d& Aj,
                            std::tuple<Eigen::Vector3d, Eigen::Vector3d>& gh) {
    // eq 77
    Eigen::Matrix3d combined = Aj*-Ai.inverse();
    Eigen::EigenSolver<Eigen::Matrix3d> eigensolver(combined);
    Eigen::Vector3cd eigenvalues = eigensolver.eigenvalues();
    std::vector<double> real_eigens;
    for(long int idx = 0; idx < eigenvalues.rows(); idx++) {
        if(eigenvalues(idx).imag() == 0) {
            real_eigens.push_back(eigenvalues(idx).real());
        }
    }
    Eigen::Vector3d g, h;
    for(const auto& e_val : real_eigens) {
        if(IntersectConics(Ai, Aj, e_val, gh)) {
            std::tie(g, h) = gh;
            return true;
        }
    }
    return false;
}

bool ChooseIntersection(const std::tuple<Eigen::Vector3d, Eigen::Vector3d>& gh, 
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
    Eigen::Vector3d gIntersect = g.cross(lineOfCenters);
    // get point where lineOfCenters and h intersect
    Eigen::Vector3d hIntersect = h.cross(lineOfCenters);

    // normalize gIntersect
    gIntersect = gIntersect.array()/gIntersect.array()[2];

    // normalize hIntersect
    hIntersect = hIntersect.array()/hIntersect.array()[2];

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
    Eigen::Matrix3d envelope = getMatrixAdjugate(locus);
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

bool computeCraterTriadInvariants(Conic& A, Conic& B, Conic& C,
                                  std::vector<double>& invariants) {
    std::tuple<Eigen::Vector3d, Eigen::Vector3d> gh_ij, gh_jk, gh_ki;
    Eigen::Vector3d lij, ljk, lki;
    double invA, invB, invC;
    Eigen::Matrix3d locusA = A.getLocus();
    Eigen::Matrix3d locusB = B.getLocus();
    Eigen::Matrix3d locusC = C.getLocus();
    Eigen::Vector2d centerA = A.getCenter();
    Eigen::Vector2d centerB = B.getCenter();
    Eigen::Vector2d centerC = C.getCenter();
    
    if(!IntersectionLines(locusA, locusB, gh_ij)) {
        std::cerr<<"IntersectionLines error ij\n";
        return false;
    }
    if(!IntersectionLines(locusB, locusC, gh_jk)) {
        std::cerr<<"IntersectionLines error jk\n";
        return false;
    }
    if(!IntersectionLines(locusC, locusA, gh_ki)) {
        std::cerr<<"IntersectionLines error ki\n";
        return false;
    }
    if(!ChooseIntersection(gh_ij, centerA, centerB, lij)) {
        std::cerr<<"ChooseIntersection error ij\n";
        return false;
    }
    if(!ChooseIntersection(gh_jk, centerB, centerC, ljk)) {
        std::cerr<<"ChooseIntersection error jk\n";
        return false;
    }
    if(!ChooseIntersection(gh_ki, centerC, centerA, lki)) {
        std::cerr<<"ChooseIntersection error ki\n";
        return false;
    }
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
    invariants.push_back(invA);
    invariants.push_back(invB);
    invariants.push_back(invC);
    return true;
}

Eigen::Matrix3d getENUFrame(const Eigen::Vector3d& surface_point) {
    Eigen::Vector3d u_north_pole(3), u_surface_point(3), ui, ei, ni;
    Eigen::Matrix3d enu_frame(3, 3);
    double eps = 1e-6;
    u_north_pole = getNorthPoleUnitVector();
    normalizeVector(surface_point, u_surface_point);
    // std::cout << "ENU: \n" << surface_point << "\n" << u_surface_point << std::endl;
    if((u_north_pole - u_surface_point).norm() < eps) {
        ui = u_north_pole;
        ei << 1, 0, 0;
        ni << 0, 1, 0;
    }
    else if((-u_north_pole - u_surface_point).norm() < eps) {
        ui = -u_north_pole;
        ei << -1, 0, 0;
        ni << 0, -1, 0;
    }
    else {
        ui = u_surface_point;
        ei = u_north_pole.cross(ui);
        ni = ui.cross(ei);
    }
    enu_frame << ei, ni, ui;
    return enu_frame;
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

Eigen::Matrix3d pointCameraInDirection(const Eigen::Vector3d& camera_position, const Eigen::Vector3d& desired_location) {
    Eigen::Matrix3d T_m2c = getENUFrame(desired_location - camera_position);
    Eigen::Matrix3d z_rot;
    eulerToQuaternion(0., 0., M_PI, z_rot);
    return z_rot * T_m2c.transpose();
}

Eigen::Quaterniond eulerToQuaternion(const double roll,
                                     const double pitch,
                                     const double yaw) {
    Eigen::AngleAxisd rollAngle(roll, Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd pitchAngle(pitch, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd yawAngle(yaw, Eigen::Vector3d::UnitZ());

    Eigen::Quaterniond q = yawAngle * pitchAngle * rollAngle;
    return q;
}

void eulerToQuaternion(const double roll,
                       const double pitch,
                       const double yaw,
                       Eigen::Matrix3d& dcm) {
    Eigen::AngleAxisd rollAngle(roll, Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd pitchAngle(pitch, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd yawAngle(yaw, Eigen::Vector3d::UnitZ());

    dcm = yawAngle * pitchAngle * rollAngle;
}
