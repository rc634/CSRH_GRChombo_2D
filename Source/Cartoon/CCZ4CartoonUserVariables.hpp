/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CCZ4VARIABLES_HPP
#define CCZ4VARIABLES_HPP

#include <algorithm>
#include <array>
#include <string>

/// This enum gives the index of the CCZ4 variables on the grid
enum
{
    c_chi,

    c_h11,
    c_h12,
    c_h22,
    c_hww,

    c_K,

    c_A11,
    c_A12,
    c_A22,
    c_Aww,

    c_Theta,

    c_Gamma1,
    c_Gamma2,

    c_lapse,

    c_shift1,
    c_shift2,

    c_B1,
    c_B2,

    c_phi,
    c_phi_Im,
    c_Pi,
    c_Pi_Im,
    NUM_CCZ4_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_CCZ4_VARS> ccz4_variable_names = {
    "chi",

    "h11",    "h12",    "h22",

    "hww",

    "K",

    "A11",    "A12",    "A22",

    "Aww",

    "Theta",

    "Gamma1", "Gamma2",

    "lapse",

    "shift1", "shift2",

    "B1",     "B2",
    
    "phi",  "phi_Im",

    "Pi",   "Pi_Im"};
} // namespace UserVariables

#endif /* CCZ4VARIABLES_HPP */
