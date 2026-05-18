/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXPANSIONSCALAR_HPP_
#define EXPANSIONSCALAR_HPP_

#include "ADMConformalVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"
#include <array>
#include <cmath>

template <class deriv_t = FourthOrderDerivatives>
class ExpansionScalar
{
  protected:
    deriv_t m_deriv;
    double m_dx;
    std::array<double, CH_SPACEDIM> m_center;

    template <class data_t>
    using Vars = CCZ4CartoonVars::VarsWithGauge<data_t>;

  public:
    ExpansionScalar(const double a_dx,
                    const std::array<double, CH_SPACEDIM> &a_center)
        : m_deriv(a_dx), m_dx(a_dx), m_center(a_center) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto vars = current_cell.template load_vars<Vars>();

        // Cylindrical radius: in the cartoon method the y-coordinate is the
        // radial distance from the symmetry axis.  Clamp to avoid 1/0 on-axis.
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const data_t radius = simd_max(coords.y, 1e-6);

        // Sweep s^i through 32 directions from [1,0] to [-1,0] (half circle),
        // with step 16 giving [0,1].  For each direction compute the full
        // physical expansion Theta = A_ss/chi + s^y/rho and keep the maximum.
        const int s_num = 32;
        data_t AssMax = 0.0;

        for (int k = 0; k < s_num; k++)
        {
            double angle = k * M_PI / s_num;
            Tensor<1, data_t> si;
            si[0] = cos(angle);
            si[1] = sin(angle);

            // Physical length squared: h_ij s^i s^j / chi
            data_t ss = 0.0;
            FOR2(i, j)
            {
                ss += vars.h[i][j] * si[i] * si[j];
            }
            ss /= vars.chi;

            // Normalise s^i to a unit vector before contracting with A_ij
            data_t norm = sqrt(ss);
            Tensor<1, data_t> si_hat;
            FOR1(i) { si_hat[i] = si[i] / norm; }

            data_t A_ss = 0.0;
            FOR2(i, j)
            {
                A_ss += vars.A[i][j] * si_hat[i] * si_hat[j];
            }

            // curved infinitessimal planes
            data_t Div_s = 0.; //something (radius);

            AssMax = simd_max(AssMax, A_ss / vars.chi + Div_s);
        }

        data_t Q = AssMax - (2.0 / 3.0) * vars.K;

        current_cell.store_vars(Q, c_Q);
    }
};

#endif /* EXPANSIONSCALAR_HPP_ */
