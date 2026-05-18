/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(EFFECTIVEPOTENTIAL_HPP_)
#error "This file should only be included through EffectivePotential.hpp"
#endif

#ifndef EFFECTIVEPOTENTIAL_IMPL_HPP_
#define EFFECTIVEPOTENTIAL_IMPL_HPP_

#include "simd.hpp"

template <class data_t>
void EffectivePotential::compute(Cell<data_t> current_cell) const
{
    // copy data from chombo gridpoint into local variables
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

    // Get the coordinates
    const Coordinates<data_t> coords(current_cell, m_deriv.m_dx, m_center);

    const data_t x = coords.x;
    const double y = coords.y;

    // floor on chi
    const data_t chi = simd_max(vars.chi, 1e-4);

    // spacetime metric
    Tensor<2, data_t> gd = 0.;
    FOR(i, j) { gd[i][j] = vars.h[i][j] / chi; }
    data_t gww = vars.hww / chi;

    data_t proper_dist;
    data_t beta_sq;
    data_t lapse_sq;

    proper_dist = sqrt(gd[1][1] * y * y - 2 * gd[1][2] * x * y + gd[2][2] * x * x);
    beta_sq = vars.shift[0] * vars.shift[0] + vars.shift[1] * vars.shift[1] + vars.shift[2] * vars.shift[2];
    lapse_sq = vars.lapse * vars.lapse;

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(proper_dist, c_proper_dist);
    current_cell.store_vars(lapse_sq, c_lapse_sq);
    current_cell.store_vars(beta_sq, c_beta_sq);
}

#endif /* EFFECTIVEPOTENTIAL_IMPL_HPP_ */