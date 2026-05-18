/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <array>

#ifndef THINSHELLPARAMS_HPP_
#define THINSHELLPARAMS_HPP_

//! A structure for the input params for the boson star
struct ThinShell_params_t
{
    std::array<double, CH_SPACEDIM> star_centre; //!< coordinates of the centre of the star
    std::string base_path;
    double BS_frequency;
};

#endif /* THINSHELLPARAMS_HPP_ */
