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
Conic::Conic(const std::tuple<float, float, float, float, float> geom_tup) {
    semimajor_axis_ = std::get<0>(geom_tup);
    semiminor_axis_ = std::get<1>(geom_tup);
    x_center_       = std::get<2>(geom_tup);
    y_center_       = std::get<3>(geom_tup);
    angle_          = std::get<4>(geom_tup);
}
Conic::Conic(const std::vector<float> geom_vec) {
    setGeometricParameters(geom_vec);
}

void Conic::setGeometricParameters(const std::vector<float> geom_vec) {
    semimajor_axis_ = geom_vec[0];
    semiminor_axis_ = geom_vec[1];
    x_center_       = geom_vec[2];
    y_center_       = geom_vec[3];
    angle_          = geom_vec[4];
}

void Conic::NormalizeImplicitParameters(std::vector<float>& impl_params) {
    auto vecNorm = vectorNorm(impl_params);
    // std::cout << "Vector norm: " << vecNorm << std::endl;
    for(auto& element : impl_params) {
        element /= vecNorm;
    }
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

std::vector<float> Conic::getGeom() {
    return {semimajor_axis_, semiminor_axis_, x_center_, y_center_, angle_};
}

void Conic::setFromImplicit(const std::vector<float> impl_params) {
    std::vector<float> geom_params = impl_to_geom(impl_params);
    setGeometricParameters(geom_params);
}

std::vector<float> Conic::impl_to_geom(const std::vector<float> impl_params){
    float A, B, C, D, E, F, B2_minus_4AC, xc, yc, phi;
    float numerator, denominator_a, denominator_b;
    float semimajor_axis, semiminor_axis, amc2, b2;
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
        // phi = (1 / 2) * np.arccot((A - C) / B)
        phi = (1 / 2) * atan2(B, A - C);
        if(A>C) {
            phi += M_PI / 2;
        }
    // return geom: 5x1 vector
    }
    return {semimajor_axis, semiminor_axis, xc, yc, phi};
}
