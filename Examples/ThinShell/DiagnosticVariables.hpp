/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_Ham,

    c_Mom1,
    c_Mom2,

    // sqrt(Mom_1^2 + Mom_2^2)
    c_Mom,

    c_Weyl4_Re,
    c_Weyl4_Im,

    c_Madm,
    c_Padm,

    c_rho,
    c_rho_ADM, // basically rho * sqrt(gamma)

    c_Sx,
    c_Sy,

    c_N, 
    c_mod_phi,

    c_dA11,
    c_dA12,
    c_dA22,
    c_dAww,
    c_dK,

    c_proper_dist,
    c_lapse_sq,
    c_beta_sq,

    c_dt_mod_phi,
    c_gamma_tt,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",

    "Mom1",     "Mom2",     "Mom",

    "Weyl4_Re", "Weyl4_Im",

    "M_adm",    "P_adm",

    "rho",  "rho_ADM",  "Sx",  "Sy",

    "N", "mod_phi",

    "dA11", "dA12", "dA22", "dAww", "dK", 

    "proper_dist", "lapse_sq", "beta_sq",

    "dt_mod_phi", "gamma_tt",
    
    };
} // namespace DiagnosticVariables

#endif /* DIAGNOSTICVARIABLES_HPP */
