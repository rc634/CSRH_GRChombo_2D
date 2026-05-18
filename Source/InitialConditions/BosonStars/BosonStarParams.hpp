/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <array>

#ifndef BOSONSTARPARAMS_HPP_
#define BOSONSTARPARAMS_HPP_

//! A structure for the input params for the boson star / binary configuration
struct BosonStar_params_t
{
    // Boson star properties
    int    BS_gridpoints;  // ODE grid points for the 1D solution
    double BS_amplitude;   // central amplitude of the scalar field
    double BS_phase;       // initial phase of the complex field
    int    BS_eigen;       // radial eigenstate (0 = ground state)

    // Binary / configuration
    bool   binary;                                 // true = binary, false = single object
    double binary_separation;                      // axial separation of the two centres
    double binary_rapidity;                        // Lorentz rapidity of each object
    std::array<double, CH_SPACEDIM> binary_centre; // centre of the binary (midpoint)
    std::array<int, 2> binary_type;                // per-component type: 0=BS, 1=BH
    double BH_mass;                                // mass of any BH component
    bool   BH_override;                            // if true, use BH_mass instead of sol.get_ADMmass() for BH pixels

    bool print_asymptotics;
    double Newtons_constant;

    // legacy / unused
    double mass_ratio;
    int n_power;
    int id_choice;
    double radius_width1;
    double radius_width2;
    int conformal_factor_power;
    bool antiboson;
};

#endif /* BOSONSTARPARAMS_HPP_ */
