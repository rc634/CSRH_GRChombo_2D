/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef NOETHERCHARGE_HPP_
#define NOETHERCHARGE_HPP_

#include "ComplexScalarField.hpp"
#include "ADMConformalVars.hpp" // needed for CCz4 and matter variables
#include "FourthOrderDerivatives.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"

//! Calculates the Noether Charge integrand values and the modulus of the
//! complex scalar field on the grid
template <class deriv_t = FourthOrderDerivatives>
class NoetherCharge
{
protected:
    deriv_t m_deriv;

    // Need matter variables and chi
    template <class data_t> using Vars
                                = CCZ4CartoonVars::VarsWithGauge<data_t>;

public:

    NoetherCharge(const double a_dx) : m_deriv(a_dx) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // load vars locally
        const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
        // const auto matter_vars = current_cell.template load_vars<MatterVars>();
        const auto advec_csf =
            this->m_deriv.template advection<Vars>(current_cell, vars.shift);

        Coordinates<data_t> coords(current_cell, this->m_deriv.m_dx);

        using namespace TensorAlgebra;

        // calculate Noether charge
        data_t N = pow(vars.chi, -1.5) * (vars.phi_Im
            * vars.Pi - vars.phi * vars.Pi_Im);

        //Calculate modulus of the scalar field
        data_t mod_phi = sqrt(vars.phi * vars.phi
                            + vars.phi_Im * vars.phi_Im);

        //Calculate time derivative phi^2
        data_t phi_Re_t = advec_csf.phi - vars.lapse * vars.Pi;
        data_t phi_Im_t = advec_csf.phi_Im - vars.lapse * vars.Pi_Im;
        data_t dt_A_sq = abs(2 * phi_Re_t * vars.phi + 2 * phi_Im_t * vars.phi_Im);

        //Calculate time derivative of g_{tt}
        data_t lapse_t = 1.0 * advec_csf.lapse - 2.0 * pow(vars.lapse, 1.0) * (vars.K - 2. * vars.Theta);
        Tensor<1, data_t> shift_t;
        FOR(i)
        {
            shift_t[i] = 0.75 * vars.B[i];
        }

        data_t abs_gamma_tt;
        Tensor<2, data_t> h_t;
        const int dI = CH_SPACEDIM - 1;
        const int nS = GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere
        const double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
        const double two_over_gr_spacedim = 2. * one_over_gr_spacedim;

        auto h_UU = compute_inverse_sym(vars.h);
        auto h_UU_ww = 1. / vars.hww;

        data_t divshift = compute_trace(d1.shift) +
                  nS * vars.shift[dI] /
                  coords.y; //!< includes the higher D contribution

        data_t trA = compute_trace(vars.A, h_UU);
        trA += nS * h_UU_ww * vars.Aww;

        Tensor<2, data_t> A_TF;
        FOR(i, j)
        {
            A_TF[i][j] = vars.A[i][j] - trA * vars.h[i][j] * one_over_gr_spacedim;
        }

        FOR(i, j)
        {
            h_t[i][j] = advec_csf.h[i][j] - 2.0 * vars.lapse * A_TF[i][j] - two_over_gr_spacedim * vars.h[i][j] * divshift;
            FOR(k)
                {
                h_t[i][j] +=
                        vars.h[k][i] * d1.shift[k][j] + vars.h[k][j] * d1.shift[k][i];
                }
        }

        data_t gamma_tt = -2.0 * vars.lapse * lapse_t;
        FOR(i, j)
        {
            gamma_tt += shift_t[i] * vars.h[i][j] * vars.shift[j] + vars.shift[i] * (h_t[i][j] * shift_t[j] + vars.h[i][j] * shift_t[j]);
        }

        abs_gamma_tt = abs(gamma_tt);

        current_cell.store_vars(N, c_N);
        current_cell.store_vars(mod_phi, c_mod_phi);
        // current_cell.store_vars(dt_A_sq, c_dt_mod_phi);
        // current_cell.store_vars(abs_gamma_tt, c_gamma_tt);
    }
};

#endif /* NOETHERCHARGE_HPP_ */
