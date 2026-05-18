/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXPHIANDCHIEXTRACTIONTAGGINGCRITERION_HPP_
#define COMPLEXPHIANDCHIEXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ComplexScalarField.hpp"
#include "SimulationParametersBase.hpp"
#include "Tensor.hpp"
#include "Coordinates.hpp"
#include "DebuggingTools.hpp"

class ComplexPhiAndChiExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const bool m_activate_extraction;
    const extraction_params_t m_params;
    const int m_level;
    const double m_threshold_phi;
    const double m_threshold_chi;

    template <class data_t>
    using MatterVars = typename ComplexScalarField<>::template Vars<data_t>;

  public:
    ComplexPhiAndChiExtractionTaggingCriterion(const double a_dx,
        const int a_level, const extraction_params_t a_params,
        const double a_threshold_phi, const double a_threshold_chi, const bool activate_extraction = false)
        : m_dx(a_dx), m_deriv(a_dx), m_params(a_params), m_level(a_level),
        m_threshold_phi(a_threshold_phi), m_threshold_chi(a_threshold_chi), m_activate_extraction(activate_extraction) {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {

        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);
        const auto d2 = m_deriv.template diff2<MatterVars>(current_cell);

        auto chi_val = current_cell.load_vars(c_chi);
        Tensor<1, data_t> d1chi;
        FOR(idir) m_deriv.diff1(d1chi, current_cell, idir, c_chi);

        data_t sum_d2_sq = 0., sum_d1_sq = 0.;
        FOR2(idir, jdir)
            sum_d2_sq += d2.phi[idir][jdir] * d2.phi[idir][jdir]
                       + d2.phi_Im[idir][jdir] * d2.phi_Im[idir][jdir]
                       + d2.Pi[idir][jdir] * d2.Pi[idir][jdir]
                       + d2.Pi_Im[idir][jdir] * d2.Pi_Im[idir][jdir];
        FOR(idir)
            sum_d1_sq += d1.phi[idir] * d1.phi[idir]
                       + d1.phi_Im[idir] * d1.phi_Im[idir]
                       + d1.Pi[idir] * d1.Pi[idir]
                       + d1.Pi_Im[idir] * d1.Pi_Im[idir];
        data_t d2_phi_ratio = sqrt(sum_d2_sq) / (sqrt(sum_d1_sq) + 1e-10);

        // Amplitude weight: suppresses phi tagging where |phi| is small.
        // The exponential tail has constant d2/d1 ratio independent of amplitude,
        // so without this the criterion fires uniformly outside the star.
        data_t phi_Re = current_cell.load_vars(c_phi);
        data_t phi_Im = current_cell.load_vars(c_phi_Im);
        data_t mod_phi_sq = phi_Re * phi_Re + phi_Im * phi_Im;
        data_t phi_weight = mod_phi_sq / (mod_phi_sq + 1e-6);

        data_t chi_grad_sq = 0.;
        FOR(idir) chi_grad_sq += d1chi[idir] * d1chi[idir];
        data_t safe_chi = simd_max(chi_val, 1e-10);

        data_t criterion_chi_puncture = 200. * m_dx * sqrt(chi_grad_sq) / (safe_chi * safe_chi * m_threshold_chi);
        data_t criterion_phi = 20. * m_dx * sqrt(d2_phi_ratio) / m_threshold_phi * phi_weight;

        data_t criterion = simd_max(criterion_phi, criterion_chi_puncture);

        if (m_activate_extraction)
        {
        for (int iradius = 0; iradius < m_params.num_extraction_radii;
             ++iradius)
        {
            // regrid if within extraction level and not at required refinement
            if (m_level < m_params.extraction_levels[iradius])
            {
                const Coordinates<data_t> coords(current_cell, m_dx,
                                                 m_params.extraction_center);
                const data_t r = coords.get_radius();
                // add a 20% buffer to extraction zone so not too near to
                // boundary
                auto regrid = simd_compare_lt(
                    r, 1.2 * m_params.extraction_radii[iradius]);
                criterion = simd_conditional(regrid, 100.0, criterion);
            }
        }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* COMPLEXPHIANDCHIEXTRACTIONTAGGINGCRITERION_HPP_ */
