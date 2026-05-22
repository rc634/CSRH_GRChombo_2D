/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTAR_HPP_)
#error "This file should only be included through BosonStar.hpp"
#endif

#ifndef BOSONSTAR_IMPL_HPP_
#define BOSONSTAR_IMPL_HPP_

#include "BosonStarSolution.hpp" //for BosonStarSolution class
#include "BSPixel.hpp"
#include <iomanip>
#include <sstream>

inline BosonStar::BosonStar(BosonStar_params_t a_params_BosonStar,
                            Potential::params_t a_params_potential,
                            double a_G_Newton, double a_dx, int a_verbosity)
    : m_dx(a_dx), m_G_Newton(a_G_Newton),
      m_params_BosonStar(a_params_BosonStar),
      m_params_potential(a_params_potential), m_verbosity(a_verbosity)
{
}

void BosonStar::print_star_info() const
{
    std::ostringstream msg;
    msg << "  Boson star:  M_ADM = " << std::fixed << std::setprecision(6)
        << m_1d_sol.get_ADMmass()
        << "   phi_c = " << m_1d_sol.get_p_interp(0.);
    std::string s = msg.str();
    const int box_width = 52;
    if ((int)s.size() < box_width)
        s.resize(box_width, ' ');
    pout() << "  +----------------------------------------------------+" << endl;
    pout() << "  |" << s << "|" << endl;
    pout() << "  +----------------------------------------------------+" << endl;
}

void BosonStar::compute_1d_solution(const double max_r)
{
    try
    {
        // set initial parameters and then run the solver (didnt put it in the
        // constructor)
        pout() << "Setting initial conditions" << endl;
        m_1d_sol.set_initialcondition_params(m_params_BosonStar,
                                             m_params_potential, max_r);
        pout() << "run m_1d_sol.main()" << endl;
        m_1d_sol.main();
        pout() << "completed m_1d_sol.main()" << endl;
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}



// Compute the value of the initial vars on the grid
template <class data_t> void BosonStar::compute(Cell<data_t> current_cell) const
{
    if (m_params_BosonStar.binary) {
        compute_BS_binary(current_cell);
    } else {
        compute_BH(current_cell);
    }

}

template <class data_t> void BosonStar::compute_BH(Cell<data_t> current_cell) const
{
    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    current_cell.load_vars(vars);

    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_BosonStar.binary_centre);

    BSPixel bh1;
    bh1.is_BH = true; // must be set before init
    bh1.m = 1.;
    bh1.init(1., coords, m_1d_sol, m_params_BosonStar);

    vars.chi      = bh1.chi;
    vars.lapse    = sqrt(vars.chi);
    vars.shift[0] = bh1.shift_x;
    vars.K        = bh1.trK;

    FOR2(i,j) {
      vars.h[i][j] = bh1.gamma[i][j] * vars.chi;
      vars.A[i][j] = (bh1.K[i][j] - (vars.K/3.) * bh1.gamma[i][j]) * vars.chi;
    }
    vars.hww = bh1.gamma[2][2] * vars.chi;
    vars.Aww = (bh1.K[2][2] - (vars.K/3.) * bh1.gamma[2][2]) * vars.chi;

    current_cell.store_vars(vars);
}

// Compute the value of the initial vars on the grid
template <class data_t> void BosonStar::compute_BS_binary(Cell<data_t> current_cell) const
{
    // Lets goooo

    CCZ4CartoonVars::VarsWithGauge<data_t> vars;
    //  MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)

    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below
    //Coordinates<data_t> coords(current_cell, m_dx,
    //                           m_params_BosonStar.star_centre);


    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_BosonStar.binary_centre);

     // Import binary configuration params
    double rapidity  = m_params_BosonStar.binary_rapidity;
    bool binary      = m_params_BosonStar.binary;
    double M         = m_params_BosonStar.BH_mass;
    double separation = m_params_BosonStar.binary_separation;
    // int conformal_power = m_params_BosonStar.conformal_factor_power;
    int initial_data_choice = m_params_BosonStar.id_choice;
    bool print_asymptotics = m_params_BosonStar.print_asymptotics;

    // * - * - * - * - * - * - * - * - * - * - * - * -

    // Here we define the full state of a single boosted BS or BH
    // for the current gridpoint or pixel

    // define two blanc star pixels
    BSPixel bs1;
    BSPixel bs2;

    // binary_type: 0 = boson star, 1 = black hole (per component)
    bs1.is_BH = (m_params_BosonStar.binary_type[0] == 1);
    bs2.is_BH = (m_params_BosonStar.binary_type[1] == 1);

    // if we have a black hole we can set the mass, To DO!

    // init must come after setting is_BH !
    bs1.init(1., coords, m_1d_sol, m_params_BosonStar);
    bs2.init(-1., coords, m_1d_sol, m_params_BosonStar);

    // * - * - * - * - * - * - * - * - * - * - * - * -


    // Here we use Thomas Helfer's trick and find the corresponding fixed values to be substracted in the initial guess
    // is kroneka delta in infinite separation
    double kroneker_helfer[3][3] = {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
    double gamma[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double Kij[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};


    // superimpose metrics and curvatures
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        gamma[i][j] = bs1.gamma[i][j]
                    + bs2.gamma[i][j]
                    - kroneker_helfer[i][j];
        Kij[i][j] = bs1.K[i][j] + bs2.K[i][j];
      }
      // DIAGONAL METRIC ASSUMED!
      vars.K += Kij[i][i]/gamma[i][i];
    }

    // GAUGE
    // take vars from pixel sum
    vars.chi = pow(gamma[0][0]*gamma[1][1]*gamma[2][2],-1./3.);
    // pre-collapsed lapse
    vars.lapse = sqrt(vars.chi);
    vars.shift[0] = bs1.shift_x + bs2.shift_x;

    // METRIC / CURVATURE
    FOR2(i,j) {
      vars.h[i][j] = gamma[i][j] * vars.chi;
      vars.A[i][j] = (Kij[i][j] - (vars.K/3.) * gamma[i][j]) * vars.chi;
    }
    vars.hww = gamma[2][2] * vars.chi;
    vars.Aww = (Kij[2][2] - (vars.K/3.) * gamma[2][2]) * vars.chi;

    // MATTER
    vars.phi = bs1.phi + bs2.phi;
    vars.phi_Im = bs1.phi_Im + bs2.phi_Im;
    vars.Pi = bs1.Pi + bs2.Pi;
    vars.Pi_Im = bs1.Pi_Im + bs2.Pi_Im;

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}



// // Compute the value of the initial vars on the grid
// template <class data_t> void BosonStar::compute_GERRYTEST(Cell<data_t> current_cell) const
// {
//     // Lets goooo

//     CCZ4CartoonVars::VarsWithGauge<data_t> vars;
//     //  MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
//     // Load variables (should be set to zero if this is a single BS)

//     current_cell.load_vars(vars);
//     // VarsTools::assign(vars, 0.); // Set only the non-zero components below
//     //Coordinates<data_t> coords(current_cell, m_dx,
//     //                           m_params_BosonStar.star_centre);


//     Coordinates<data_t> coords(current_cell, m_dx,
//                                m_params_BosonStar.binary_centre);

//      // Import binary configuration params
//     double rapidity  = m_params_BosonStar.binary_rapidity;
//     bool binary      = m_params_BosonStar.binary;
//     double M         = m_params_BosonStar.BH_mass;
//     double separation = m_params_BosonStar.binary_separation;
//     // int conformal_power = m_params_BosonStar.conformal_factor_power;
//     int initial_data_choice = m_params_BosonStar.id_choice;
//     bool print_asymptotics = m_params_BosonStar.print_asymptotics;

//     // * - * - * - * - * - * - * - * - * - * - * - * -

//     // hard coded test metric 
//     double d = 5.;
//     double m = 0.05;
//     double r1 = sqrt((coords.x-d/2.)*(coords.x-d/2.) + coords.y*coords.y);
//     double r2 = sqrt((coords.x+d/2.)*(coords.x+d/2.) + coords.y*coords.y);
//     double psi = 1. + 0.5 * m * ( 1./r1 + 1./r2 );

//     // * - * - * - * - * - * - * - * - * - * - * - * -


//     // Here we use Thomas Helfer's trick and find the corresponding fixed values to be substracted in the initial guess
//     // is kroneka delta in infinite separation
//     double kroneker_helfer[3][3] = {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
//     double gamma[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
//     double Kij[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};


//     // superimpose metrics and curvatures
//     for (size_t i = 0; i < 3; i++) {
//       for (size_t j = 0; j < 3; j++) {
//         gamma[i][j] = pow(psi,4) * kroneker_helfer[i][j];
//         Kij[i][j] = -pow(psi,4) * kroneker_helfer[i][j];
//       }
//     }
//     // DIAGONAL METRIC ASSUMED!
//     vars.K = - 3. * pow(psi,4);

//     // GAUGE
//     // take vars from pixel sum
//     vars.chi = pow(gamma[0][0]*gamma[1][1]*gamma[2][2],-1./3.);
//     // pre-collapsed lapse
//     vars.lapse = sqrt(vars.chi);

//     // METRIC / CURVATURE
//     FOR2(i,j) {
//       vars.h[i][j] = gamma[i][j] * vars.chi;
//       vars.A[i][j] = (Kij[i][j] - (vars.K/3.) * gamma[i][j]) * vars.chi;
//     }
//     vars.hww = gamma[2][2] * vars.chi;
//     vars.Aww = (Kij[2][2] - (vars.K/3.) * gamma[2][2]) * vars.chi;

//     // no matter

//     // Store the initial values of the variables
//     current_cell.store_vars(vars);
// }


#endif /* BOSONSTAR_IMPL_HPP_ */
