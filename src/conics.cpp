#include <math.h>

#include "conics.h"

Conic::Conic(const float semimajor_axis, 
             const float semiminor_axis, 
             const float x_center, 
             const float y_center, 
             const float angle) {
    semimajor_axis_ = semimajor_axis;
    semiminor_axis_ = semiminor_axis;
    x_center_ = x_center;
    y_center_ = y_center;
    angle_ = angle;
}
Conic::Conic(const std::tuple<float, float, float, float, float>& geom_tup) {
    semimajor_axis_ = std::get<0>(geom_tup);
    semiminor_axis_ = std::get<1>(geom_tup);
    x_center_       = std::get<2>(geom_tup);
    y_center_       = std::get<3>(geom_tup);
    angle_          = std::get<4>(geom_tup);
}
Conic::Conic(const std::vector<float>& geom_vec) {
    setGeometricParameters(geom_vec);
}

void Conic::setGeometricParameters(const std::vector<float>& geom_vec) {
    semimajor_axis_ = geom_vec[0];
    semiminor_axis_ = geom_vec[1];
    x_center_       = geom_vec[2];
    y_center_       = geom_vec[3];
    angle_          = geom_vec[4];
}

void Conic::setImplicitParameters(const std::vector<float>& impl_params) {
    std::vector<float> geom_params = Implicit2Geom(impl_params);
    setGeometricParameters(geom_params);
}

void Conic::setLocus(const Eigen::Matrix3f& locus) {
    std::vector<float> geom_params = Locus2Geom(locus);
    setGeometricParameters(geom_params);
}

void Conic::NormalizeImplicitParameters(std::vector<float>& impl_params) {
    auto vecNorm = vectorNorm(impl_params);
    // std::cout << "Vector norm: " << vecNorm << std::endl;
    for(auto& element : impl_params) {
        element /= vecNorm;
    }
}

std::vector<float> Conic::getGeom() {
    return {semimajor_axis_, semiminor_axis_, x_center_, y_center_, angle_};
}

std::vector<float> Conic::getImplicit() {
    return Geom2Implicit();
}

Eigen::Matrix3f Conic::getLocus() {
    return Geom2Locus();
}

Eigen::Matrix3f Conic::getEnvelope() {
    return getMatrixAdjugate(Geom2Locus());
}

std::vector<float> Conic::Locus2Implicit(const Eigen::Matrix3f& locus) {
    std::vector<float> coeff;
    coeff.reserve(6);
    coeff.push_back(locus.coeff(0,0));
    coeff.push_back(2*locus.coeff(0,1));
    coeff.push_back(locus.coeff(1,1));
    coeff.push_back(2*locus.coeff(0,2));
    coeff.push_back(2*locus.coeff(1,2));
    coeff.push_back(locus.coeff(2,2));
    return coeff;
}

std::vector<float> Conic::Implicit2Geom(const std::vector<float>& impl_params){
    float numerator, denominator_a, denominator_b;
    float B2_minus_4AC, xc, yc, phi;
    float semimajor_axis, semiminor_axis, amc2, b2;
    float A, B, C, D, E, F;
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

Eigen::Matrix3f Conic::Geom2Locus() {
    const std::vector<float> impl_params = Geom2Implicit();
    // return normalize_determinant(Implicit2Locus(impl_params))
    return Implicit2Locus(impl_params);
 }

Eigen::Matrix3f Conic::Implicit2Locus(const std::vector<float>& impl_params) {
    Eigen::Matrix3f locus_mtx(3,3);
    float A, B, C, D, E, F;
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

 std::vector<float> Conic::Locus2Geom(const Eigen::Matrix3f& locus) {
    const std::vector<float> impl_params = Locus2Implicit(locus);
    return Implicit2Geom(impl_params);
 }


std::vector<float> Conic::Geom2Implicit() {
    // unpack geometric parameters
    float a   = semimajor_axis_;
    float b   = semiminor_axis_;
    float xc  = x_center_;
    float yc  = y_center_;
    float phi = angle_;

    // perform some computations beforehand
    float a2 = pow(a, 2);
    float b2 = pow(b, 2);
    float sin_phi = sin(phi);
    float cos_phi = cos(phi);
    float xc_2 = pow(xc, 2);
    float yc_2 = pow(yc, 2);
    float sin_phi_2 = pow(sin_phi, 2);
    float cos_phi_2 = pow(cos_phi, 2);

    // Populuate each coefficient
    std::vector<float> coeff;
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

bool Conic::ConicIntersectionLines(const Eigen::Matrix3f& Aj, 
                                   std::tuple<Eigen::Vector3f, Eigen::Vector3f>& gh) {
    Eigen::Matrix3f Ai = getLocus();
    return IntersectionLines(Ai, Aj, gh);
}

bool Conic::ConicIntersectionLines(Conic& conicB, 
                                   std::tuple<Eigen::Vector3f, Eigen::Vector3f>& gh) {
    Eigen::Matrix3f Ai = getLocus();
    Eigen::Matrix3f Aj = getLocus();
    return IntersectionLines(Ai, Aj, gh);
}

Eigen::Vector2f Conic::getCenter() {
    Eigen::Vector2f center;
    center << x_center_, y_center_;
    return center;
}

Eigen::Vector2f Conic::getSemiAxes() {
    Eigen::Vector2f axes;
    axes << semimajor_axis_, semiminor_axis_;
    return axes;
}

std::vector<float> convertEigenVectorToVector(const Eigen::Vector3f& eig) {
    std::vector<float> regVec;
    regVec.resize(eig.size());
    Eigen::Vector3f::Map(&regVec[0], eig.size()) = eig;
    return regVec;
}

template <typename T>
bool vectorContainsNaN(const std::vector<T> vec) {
    return std::any_of(vec.begin(), vec.end(), [](float i){return std::isnan(i);});
}

bool vectorContainsNaN(const Eigen::Vector3f& eV) {
    std::vector<float> vec = convertEigenVectorToVector(eV);
    return vectorContainsNaN(vec);
}

Eigen::Matrix3f getMatrixAdjugate(const Eigen::Matrix3f& mtx) {

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
    Eigen::Matrix3f mtx_adj = Eigen::Matrix3f(3, 3);
    mtx_adj << mtx_adj11, mtx_adj12, mtx_adj13,
               mtx_adj21, mtx_adj22, mtx_adj23, 
               mtx_adj31, mtx_adj32, mtx_adj33;
    return mtx_adj;
}

bool IntersectConics(const Eigen::Matrix3f& Ai, 
                     const Eigen::Matrix3f& Aj, 
                     const float eig,
                     std::tuple<Eigen::Vector3f, Eigen::Vector3f>& gh) {
    float eps = 1e-16;
    Eigen::Matrix3f Bij(3, 3), Bij_star(3, 3);
    Eigen::Vector3f bkk_eig(3);
    Bij = eig*Ai + Aj;
    Bij_star = getMatrixAdjugate(Bij);

    bkk_eig = Bij_star.diagonal();
    std::vector<float> bkk = convertEigenVectorToVector(bkk_eig);
    // eq 86: all diagonal entries of Bij_star are negative
    if ( std::any_of(bkk.begin(), bkk.end(), [&eps](float i){return i>eps;}) ) {
        std::cout << "bkk contains positive numbers: " << bkk_eig << std::endl;
        return false;
    }
    int min_idx = arg_min(bkk);

    // eq 87
    Eigen::Vector3f z = -Bij_star.row(min_idx)/sqrt(-bkk[min_idx]);
    Eigen::Matrix3f D = Bij + crossMatrix(z);
    Eigen::Matrix3f Dabs = D.cwiseAbs();
    
    Eigen::Matrix3f::Index maxRow,maxCol;
    Dabs.maxCoeff(&maxRow, &maxCol);

    // pg 44 of Christian, Derksen, Watkins [2020]
    // g is non-zero column of D, h is non-zero row of D
    Eigen::Vector3f g = D.col(maxCol);
    Eigen::Vector3f h = D.row(maxRow);
    
    std::vector<float> gVec = convertEigenVectorToVector(g);
    std::vector<float> hVec = convertEigenVectorToVector(h);
    if ( vectorContainsNaN(gVec) ) {
        std::cout << "g vector contains nan's: " << g << std::endl;
        return false;
    }
    if ( vectorContainsNaN(hVec) ) {
        std::cout << "hs vector contains nan's: " << h << std::endl;
        return false;
    }
    gh = {g, h};
    return true;
}

bool IntersectionLines(const Eigen::Matrix3f& Ai, 
                            const Eigen::Matrix3f& Aj,
                            std::tuple<Eigen::Vector3f, Eigen::Vector3f>& gh) {
    // eq 77
    Eigen::Matrix3f combined = Aj*-Ai.inverse();
    Eigen::EigenSolver<Eigen::Matrix3f> eigensolver(combined);
    // Eigen::Matrix3cf eigenvectors = eigensolver.eigenvectors();
    Eigen::Vector3cf eigenvalues = eigensolver.eigenvalues();
    std::vector<float> real_eigens;
    for(long int idx = 0; idx < eigenvalues.rows(); idx++) {
        if(eigenvalues(idx).imag() == 0) {
            real_eigens.push_back(eigenvalues(idx).real());
        }
    }
    printVector(real_eigens, "Real eigenvalues: ");

    Eigen::Vector3f g, h;
    // std::tuple<Eigen::Vector3f, Eigen::Vector3f> gh;
    for(const auto& e_val : real_eigens) {
        if(IntersectConics(Ai, Aj, e_val, gh)) {
            std::tie(g, h) = gh;
            std::cout << " g: " << g << std::endl;
            std::cout << " h: " << h << std::endl;
            return true;
        }
    }
    return false;
}

bool ChooseIntersection(const std::tuple<Eigen::Vector3f, Eigen::Vector3f>& gh, 
                                   const Eigen::Vector2f& centerA, 
                                   const Eigen::Vector2f& centerB,
                                   Eigen::Vector3f& l) {
    Eigen::Vector3f g, h;
    std::tie(g, h) = gh;
    // convert centers to homogeneous coordinates
    Eigen::Vector3f centerAHom = centerA.homogeneous();
    Eigen::Vector3f centerBHom = centerB.homogeneous();
    // get line connecting the two centers
    Eigen::Vector3f lineOfCenters = centerAHom.cross(centerBHom);
    // get point where lineOfCenters and g intersect
    Eigen::Vector3f gIntersect = g.cross(lineOfCenters);
    // get point where lineOfCenters and h intersect
    Eigen::Vector3f hIntersect = h.cross(lineOfCenters);

    // std::cout << " Gint: " << std::endl << gIntersect << std::endl;
    // std::cout << " Hint: " << std::endl << hIntersect << std::endl;

    // normalize gIntersect
    gIntersect = gIntersect.array()/gIntersect.array()[2];

    // normalize hIntersect
    hIntersect = hIntersect.array()/hIntersect.array()[2];

    float xmax, xmin, ymax, ymin;
    xmax = (centerA(0)>centerB(0)) ? centerA(0) : centerB(0);
    xmin = (centerA(0)<centerB(0)) ? centerA(0) : centerB(0);
    ymax = (centerA(1)>centerB(1)) ? centerA(1) : centerB(1);
    ymin = (centerA(1)<centerB(1)) ? centerA(1) : centerB(1);
    bool g_fits, h_fits;

    std::cout << "Center A:\n" << centerA << std::endl;
    std::cout << "Center B:\n" << centerB << std::endl;
    std::cout << "xvals:\n" << xmin << ", " << xmax << std::endl;
    std::cout << "yvals:\n" << ymin << ", " << ymax << std::endl;


    std::cout << " gIintersect: " << std::endl << gIntersect << std::endl;
    std::cout << " hIintersect: " << std::endl << hIntersect << std::endl;
    g_fits = gIntersect(0)>xmin && gIntersect(0)<xmax && gIntersect(1)>ymin && gIntersect(1)<ymax;
    h_fits = hIntersect(0)>xmin && hIntersect(0)<xmax && hIntersect(1)>ymin && hIntersect(1)<ymax;
    if (!(g_fits ^ h_fits)) {
        std::cerr << "G|H -> " << g_fits << "|" << h_fits << std::endl;
        std::cerr << "Both or neither g_fits and h_fits are valid." << std::endl;
        return false;
    }
    l = g_fits ? g : h;
    std::cout << "ChooseIntersection chose " << (g_fits ? "g\n" : "h\n") << l << std::endl;
    return true;
}


// def conic_intersection(locus_a, locus_b, center_a, center_b):
//     gab, hab = conic_intersection_lines(locus_a, locus_b)
//     lab = choose_conic_intersection_line(gab, hab, center_a, center_b)
//     return lab

// def compute_invariant(line1, line2, locus_i):
bool computeInvariant(const Eigen::Vector3f& line1, 
                      const Eigen::Vector3f& line2, 
                      const Eigen::Matrix3f& locus,
                      float& invariant) {
    // if(vectorContainsNaN(line1) || vectorContainsNaN(line2)) {
    //     std::cerr << "One of the lines contains a NaN." << std::endl;
    //     return false;
    // }
    Eigen::Matrix3f envelope = getMatrixAdjugate(locus);
    // the numerator can be negative
    Eigen::RowVector3f line1T = line1.transpose();
    Eigen::RowVector3f line2T = line2.transpose();
    std::cout << "Line1T: \n" << line1T << std::endl;
    Eigen::MatrixXf numerator = line1T * envelope * line2;
    float numeratord = numerator.value();
    std::cout << "Numerator: " << numeratord << std::endl;
    // the denominator is quadratic in both terms, so never negative
    Eigen::MatrixXf l1tel1 = line1T*envelope*line1;
    Eigen::MatrixXf l2tel2 = line2T*envelope*line2;
    float l1tel1d = l1tel1.value();
    float l2tel2f = l2tel2.value();
    float denominator = sqrt(l1tel1d*l2tel2f);
    // calculate the invariant, eqns 95-97
    invariant = acosh(abs(numeratord)/denominator);
    return false;
}


// def compute_invariant_triad(locus_a, locus_b, locus_c):
//     conic_a = geo.Conic(locus=locus_a)
//     conic_b = geo.Conic(locus=locus_b)
//     conic_c = geo.Conic(locus=locus_c)

//     center_a = conic_a.get_geom_array()[2:4]
//     center_b = conic_b.get_geom_array()[2:4]
//     center_c = conic_c.get_geom_array()[2:4]
//     lab = conic_intersection(locus_a, locus_b, center_a, center_b)
//     lbc = conic_intersection(locus_b, locus_c, center_b, center_c)
//     lca = conic_intersection(locus_c, locus_a, center_c, center_a)

//     // Compute invariants
//     Ja = compute_invariant(lab, lca, locus_a)
//     Jb = compute_invariant(lbc, lab, locus_b)
//     Jc = compute_invariant(lca, lbc, locus_c)
//     return np.array([Ja, Jb, Jc])