/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICAIJ_HPP_
#define DIAGNOSTICAIJ_HPP_

#include "CCZ4CartoonGeometry.hpp"
#include "CCZ4CartoonVars.hpp"
#include "CCZ4Geometry.hpp"
#include "CCZ4RHS.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

#include <array>

template <class gauge_t = MovingPunctureGauge,
          class deriv_t = FourthOrderDerivatives,
          class potential_t = Potential>
class DiagnosticAij : public CCZ4RHS<gauge_t>
{
  public:
    using CCZ4 = CCZ4RHS<gauge_t, deriv_t>;

    using params_t = typename CCZ4RHS<gauge_t, deriv_t>::params_t;

    /// CCZ4 cartoon variables
    template <class data_t> using Vars = CCZ4CartoonVars::VarsWithGauge<data_t>;

    /// CCZ4 cartoon variables
    template <class data_t>
    using Diff2Vars = CCZ4CartoonVars::Diff2VarsWithGauge<data_t>;

    /// Constructor
    DiagnosticAij<gauge_t, deriv_t, potential_t>(params_t params, double dx, double sigma, potential_t a_potential, double a_G_Newton, 
                  int formulation, double cosmological_constant=0) : 
                  CCZ4RHS<gauge_t, deriv_t>(params, dx, sigma, formulation,
                  cosmological_constant), m_potential(a_potential), m_G_Newton(a_G_Newton) 
                  {}

    // Compute function
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        using namespace TensorAlgebra;

    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec = this->m_deriv.template advection<Vars>(current_cell, vars.shift);

    Coordinates<data_t> coords(current_cell, this->m_deriv.m_dx);
    const double cartoon_coord = coords.y;

    const int dI = CH_SPACEDIM - 1;
    const int nS = GR_SPACEDIM - CH_SPACEDIM; //!< Dimensions of the transverse sphere
    const double one_over_gr_spacedim = 1. / ((double)GR_SPACEDIM);
    const double two_over_gr_spacedim = 2. * one_over_gr_spacedim;

    auto h_UU = compute_inverse_sym(vars.h);
    auto h_UU_ww = 1. / vars.hww;
    auto chris = compute_christoffel(d1.h, h_UU);
    const double one_over_cartoon_coord = 1. / cartoon_coord;
    const double one_over_cartoon_coord2 = 1. / (cartoon_coord * cartoon_coord);

    Tensor<1, data_t> chris_ww;
    FOR(i)
        {
        chris_ww[i] = one_over_cartoon_coord * (delta(i, dI) - h_UU[i][dI] * vars.hww);
        FOR(j) chris_ww[i] -= 0.5 * h_UU[i][j] * d1.hww[j];
        }
    Tensor<1, data_t> chris_contracted; //!< includes the higher D contributions
    FOR(i)
        {
            chris_contracted[i] = chris.contracted[i] + nS * h_UU_ww * chris_ww[i];
        }

    Tensor<1, data_t> Z_over_chi;
    Tensor<1, data_t> Z;
    if (this->m_formulation == CCZ4::USE_BSSN)
    {
        FOR(i) Z_over_chi[i] = 0.0;
    }
    else
    {
        FOR(i) Z_over_chi[i] = 0.5 * (vars.Gamma[i] - chris_contracted[i]);
    }
    FOR(i) Z[i] = vars.chi * Z_over_chi[i];

    auto ricci = CCZ4CartoonGeometry::compute_ricci_Z(
        vars, d1, d2, h_UU, h_UU_ww, chris, Z_over_chi, cartoon_coord);

    data_t divshift = compute_trace(d1.shift) +
                  nS * vars.shift[dI] /
                  cartoon_coord; //!< includes the higher D contribution
    data_t divshift_w = compute_trace(d1.shift) -
                    CH_SPACEDIM * vars.shift[dI] /
                    cartoon_coord; //!< combination that appears in the

    data_t Z_dot_d1lapse = compute_dot_product(Z, d1.lapse);
    data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);
    data_t dlapse_dot_dhww = compute_dot_product(d1.lapse, d1.hww, h_UU);

    data_t hUU_dlapse_cartoon = 0;
    FOR(i)
        {
            hUU_dlapse_cartoon += one_over_cartoon_coord * h_UU[dI][i] * d1.lapse[i];
        }

    Tensor<2, data_t> covdtilde2lapse;
    Tensor<2, data_t> covd2lapse; // NOTE: we compute chi * D_i D_j lapse
    FOR(k, l)
    {
    covdtilde2lapse[k][l] = d2.lapse[k][l];
    FOR(m) { covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m]; }
        covd2lapse[k][l] =
        vars.chi * covdtilde2lapse[k][l] +
        0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
        vars.h[k][l] * dlapse_dot_dchi);
    }
    data_t covd2lapse_ww = vars.chi * (0.5 * dlapse_dot_dhww + vars.hww * hUU_dlapse_cartoon) - 0.5 * vars.hww * dlapse_dot_dchi;

    data_t tr_covd2lapse = compute_trace(covd2lapse, h_UU);
    tr_covd2lapse += nS * h_UU_ww * covd2lapse_ww;

    data_t trA = compute_trace(vars.A, h_UU);
    trA += nS * h_UU_ww * vars.Aww;

    Tensor<2, data_t> A_TF;
    FOR(i, j)
    {
        A_TF[i][j] = vars.A[i][j] - trA * vars.h[i][j] * one_over_gr_spacedim;
    }
    data_t Aww_TF;
    Aww_TF = vars.Aww - trA * vars.hww * one_over_gr_spacedim;

    Tensor<2, data_t> A_UU = raise_all(A_TF, h_UU);
    data_t A_UU_ww = h_UU_ww * h_UU_ww * Aww_TF;
    data_t tr_A2 = compute_trace(A_TF, A_UU) + nS * A_UU_ww * Aww_TF;

    Tensor<2, data_t> Adot_TF;
    FOR(i, j)
    {
        Adot_TF[i][j] = -covd2lapse[i][j] + vars.chi * vars.lapse * ricci.LL[i][j];
    }
    data_t Adot_TF_ww = -covd2lapse_ww + vars.chi * vars.lapse * ricci.LLww;

    data_t trAdot = compute_trace(Adot_TF, h_UU) + nS * h_UU_ww * Adot_TF_ww;

    FOR(i, j) { Adot_TF[i][j] -= one_over_gr_spacedim * vars.h[i][j] * trAdot; }
    Adot_TF_ww -= one_over_gr_spacedim * vars.hww * trAdot;

    Tensor<2, data_t> Aij;
    data_t Aww;
    data_t K_scalar;

    FOR(i, j)
    {
        Aij[i][j] = advec.A[i][j] + Adot_TF[i][j] +
        vars.A[i][j] * (vars.lapse * (vars.K - 2 * vars.Theta) - two_over_gr_spacedim * divshift);
    FOR(k)
        {
            Aij[i][j] += A_TF[k][i] * d1.shift[k][j] + A_TF[k][j] * d1.shift[k][i];
            FOR(l)
            {
                Aij[i][j] -= 2 * vars.lapse * h_UU[k][l] * A_TF[i][k] * A_TF[l][j];
            }
        }
    }

    Aww = advec.Aww + Adot_TF_ww +
          Aww_TF * (vars.lapse *
                    ((vars.K - 2 * vars.Theta) - 2 * h_UU_ww * Aww_TF) -
                    two_over_gr_spacedim * divshift_w);

#ifdef COVARIANTZ4
data_t kappa1_lapse = this->m_params.kappa1;
#else
data_t kappa1_lapse = this->m_params.kappa1 * vars.lapse;
#endif

    K_scalar =
    advec.K +
    vars.lapse * (ricci.scalar + vars.K * (vars.K - 2 * vars.Theta)) -
    kappa1_lapse * GR_SPACEDIM * (1 + this->m_params.kappa2) * vars.Theta -
    tr_covd2lapse;

    K_scalar += -2 * vars.lapse * GR_SPACEDIM / ((double)GR_SPACEDIM - 1.) * this->m_cosmological_constant;

    // Matter contributions

    // Need to write a function that returns potential at some point.
    data_t V_of_phi = 0.0;  // note that! here is V_of_modulus_phi_squared actually. but i am lazy to modify
    data_t dVdmodulus_phi_squared = 0.0;

    m_potential.compute_potential(V_of_phi, dVdmodulus_phi_squared, vars);

    emtensorCartoon_t<data_t> emtensor = compute_SF_EM_tensor(vars, d1, h_UU, h_UU_ww, chris, nS, V_of_phi);

    Tensor<2, data_t> Sij_TF;
    data_t Sww_TF;

    FOR(i, j)
    {
        // Need to divide by chi here to compensate for factor of chi in
        // emtensor.S
        Sij_TF[i][j] = emtensor.Sij[i][j] - one_over_gr_spacedim * emtensor.S * vars.h[i][j] / vars.chi;
    }

    Sww_TF = emtensor.Sww - one_over_gr_spacedim * emtensor.S * vars.hww / vars.chi;

    FOR(i, j)
    {
        Aij[i][j] += -8.0 * M_PI * m_G_Newton * vars.chi * vars.lapse * Sij_TF[i][j];
    }

    Aww += -8.0 * M_PI * m_G_Newton * vars.chi * vars.lapse * Sww_TF;

    K_scalar += 4.0 * M_PI * m_G_Newton * vars.lapse * (emtensor.S - 3 * emtensor.rho);

    current_cell.store_vars(Aij[0][0], c_dA11);
    current_cell.store_vars(Aij[1][1], c_dA22);
    current_cell.store_vars(Aij[0][1], c_dA12);
    current_cell.store_vars(Aww, c_dAww);
    current_cell.store_vars(K_scalar, c_dK);
}

  protected:
    const double m_G_Newton;
    const potential_t m_potential;
};

#endif /* DIAGNOSTICAIJ_HPP_ */
