/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <array>

#ifndef ROTATINGBOSONSTARPARAMS_HPP_
#define ROTATINGBOSONSTARPARAMS_HPP_

//! A structure for the input params for the boson star
struct RotatingBosonStar_params_t
{
    std::array<double, CH_SPACEDIM> star_centre; //!< coordinates of the centre of the star
    std::string base_path;
    double BS_frequency;
};

#endif /* ROTATINGBOSONSTARPARAMS_HPP_ */
