// class Conic():

// void geom_to_impl(self, geom):
//     # unpack geometric parameters
//     a, b, xc, yc, phi = geom

//     # perform some computations beforehand
//     a2 = a**2
//     b2 = b**2
//     sin_phi = np.sin(phi)
//     cos_phi = np.cos(phi)

//     # Populuate each coefficient
//     coeff = np.empty((6,))
//     coeff[0] =  a2*sin_phi**2 + b2*cos_phi**2;
//     coeff[1] =  2*(b2-a2)*sin_phi*cos_phi;
//     coeff[2] =  a2*cos_phi**2 + b2*sin_phi**2;
//     coeff[3] = -2*coeff[0]*xc - coeff[1]*yc;
//     coeff[4] = -coeff[1]*xc - 2*coeff[2]*yc;
//     coeff[5] =  coeff[0]*xc**2 + coeff[1]*xc*yc + coeff[2]*yc**2 - a2*b2;

//     # normalize coefficients
//     return  self.normalize_impl_coeff(coeff)