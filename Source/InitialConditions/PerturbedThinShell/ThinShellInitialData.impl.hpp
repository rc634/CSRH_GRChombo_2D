/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(THINSHELLINITIALDATA_HPP_)
#error "This file should only be included through ThinShellInitialData.hpp"
#endif

#ifndef THINSHELLINITIALDATA_IMPL_HPP_
#define THINSHELLINITIALDATA_IMPL_HPP_

#include <cmath>
#include "ThinShellInitialData.hpp"
#include "DebuggingTools.hpp"

inline ThinShellInitialData::ThinShellInitialData(ThinShell_params_t a_params_ThinShell, double a_dx)
    : m_params_ThinShell(a_params_ThinShell), m_dx(a_dx)
{
}

void ThinShellInitialData::compute_1d_solution()
{   
    try
    {  
        thin_shell_sol.main(m_params_ThinShell.base_path);
        pout() << "Wooo I have read the spinning BS initial data!" << endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
}

// Compute the value of the initial vars on the grid
template <class data_t>
void ThinShellInitialData::compute(Cell<data_t> current_cell) const
{   
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)
    
    current_cell.load_vars(vars);
    //VarsTools::assign(vars, 0.); // Set only the non-zero components below
    
    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
        m_params_ThinShell.star_centre);

    // Star positioning
    double t = 0;
    double x = coords.x;
    double z = 0;
    double y = coords.y;
    double r = sqrt(x * x + y * y + z * z);

    double ff = m_params_ThinShell.BS_frequency;

    int m = thin_shell_sol.m;
    
    double A_val = thin_shell_sol.get_amp_interp(r, m);
    double lapse_val = thin_shell_sol.get_lapse_interp(r, m);
    double Pi_val = 2.0 * thin_shell_sol.get_pi_interp(r, m);

    vars.lapse = lapse_val;
    double beta_x = 0;
    double beta_y = 0;
    double beta_z = 0;

    vars.shift[0] += beta_x;

    double phase_ = ff * t;

    //Add on to evolution equations
    vars.phi += A_val;
    vars.phi_Im += 0.0;
    vars.Pi += 0.0;
    vars.Pi_Im += Pi_val;

    vars.chi = thin_shell_sol.get_conformal_factor_interp(r, m);

    // Initialise conformal metric
    double h[2][2] = {{1., 0.}, {0., 1.}};

    // Define initial lapse
    FOR2(i, j) vars.h[i][j] = h[i][j];
    vars.hww = 1.0;

    // Define initial trace of K and A_ij
    vars.K = 0.0;

    FOR2(i,j) vars.A[i][j] = 0.0;

    vars.Aww = 0.0;

    current_cell.store_vars(vars);
}

#endif /* ROTATINGBOSONSTAR_IMPL_HPP_ */
