/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(ROTATINGBOSONSTAR_HPP_)
#error "This file should only be included through RotatingBosonStar.hpp"
#endif

#ifndef ROTATINGBOSONSTAR_IMPL_HPP_
#define ROTATINGBOSONSTAR_IMPL_HPP_

#include <cmath>
#include "RotatingBosonStarSolution.hpp"
#include "DebuggingTools.hpp"

inline RotatingBosonStar::RotatingBosonStar(RotatingBosonStar_params_t a_params_RotatingBosonStar, double a_dx)
    : m_params_RotatingBosonStar(a_params_RotatingBosonStar), m_dx(a_dx)
{
}

void RotatingBosonStar::compute_1d_rotating_solution()
{   
    try
    {  
        rotating_BS_sol.main(m_params_RotatingBosonStar.base_path);
        pout() << "Wooo I have read the spinning BS initial data!" << endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
}

// Compute the value of the initial vars on the grid
template <class data_t>
void RotatingBosonStar::compute(Cell<data_t> current_cell) const
{   
    double theta, phi;

    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)
    
    current_cell.load_vars(vars);
    //VarsTools::assign(vars, 0.); // Set only the non-zero components below
    
    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
        m_params_RotatingBosonStar.star_centre);

    // Star positioning
    double x = coords.x;
    double y = coords.y;
    double r = sqrt(x * x + y * y); 
    theta = acos(x/r);

    // Compactified coordinate
    double xvar = r / (1. + r);

    double ff = m_params_RotatingBosonStar.BS_frequency;

    int n = rotating_BS_sol.n;
    int m = rotating_BS_sol.m;
    
    double A_val = rotating_BS_sol.get_amp_interp(xvar, theta, n, m)*sqrt(2);
    double f_val = rotating_BS_sol.get_f_interp(xvar, theta, n, m);
    double g_val = rotating_BS_sol.get_g_interp(xvar, theta, n, m);
    double l_val = rotating_BS_sol.get_l_interp(xvar, theta, n, m);
    double omega_val = rotating_BS_sol.get_omega_interp(xvar, theta, n, m);
    double dthomega_val = rotating_BS_sol.get_dthomega_interp(xvar, theta, n, m);
    double dromega_val = rotating_BS_sol.get_dromega_interp(xvar, theta, n, m);

    double lapse = sqrt(fabs(f_val));
    double beta_x = 0.0;
    double beta_y = 0.0;
    double beta_z = -sin(theta)*omega_val;

    vars.shift[0] += beta_x;
    vars.shift[1] += beta_y;
    // vars.shift[2] += beta_z;

    double g_xx_1 = (g_val * l_val) / f_val;
    double g_yy_1 = (g_val * l_val) / f_val;
    double g_zz_1 = l_val / f_val;
    double g_xy_1 = 0.0;

    //Add on to evolution equations
    vars.phi += A_val;
    vars.phi_Im += 0.0;
    vars.Pi += 0.0;
    vars.Pi_Im += -(A_val / lapse) * (ff + omega_val/r);

    //Initialise extrinsic curvature and metric with upper indices
    double KLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double K;

    // Fill them in
    gammaLL[0][0] = g_xx_1;
    gammaLL[1][1] = g_yy_1;
    gammaLL[2][2] = g_zz_1;
    // gammaLL[0][1] = g_xy_1;
    // gammaLL[1][0] = gammaLL[0][1];

    KLL[1][1] = 0.0;
    KLL[2][2] = 0.0;
    KLL[0][0] = 0.0;
    KLL[1][2] = (l_val*(-cos(theta) * sin(theta) * dthomega_val+sin(theta)*sin(theta)*(omega_val - r*dromega_val)))/(2*f_val*r*lapse);
    KLL[2][1] = KLL[1][2];
    KLL[0][2] = (l_val*(sin(theta)*sin(theta)*dthomega_val + cos(theta)*sin(theta)*(omega_val - r*dromega_val)))/(2.*f_val*r*lapse);
    KLL[2][0] = KLL[0][2];

    double chi_arg = (g_val * g_val * l_val * l_val * l_val  / (pow(f_val, 3)));
    vars.chi = pow(chi_arg, -1. / 3.);

    // Define initial lapse
    vars.lapse += lapse;

    // Define initial trace of K and A_ij
    double one_third = 1./3.;
    FOR2(i,j) vars.h[i][j] = vars.chi * gammaLL[i][j];

    vars.hww = vars.chi * g_zz_1;

    // FOR2(i,j) vars.K += KLL[i][j] * gammaUU[i][j];
    vars.K = 0.0;
    FOR2(i,j) vars.A[i][j] = vars.chi * KLL[i][j];

    vars.Aww = vars.chi * KLL[2][2];

    current_cell.store_vars(vars);
}

#endif /* ROTATINGBOSONSTAR_IMPL_HPP_ */
