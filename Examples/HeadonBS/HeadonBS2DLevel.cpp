/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */
#include <fstream>
#include <mutex>
#include <memory>

#include "HeadonBS2DLevel.hpp"

// General headers
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"

// For tag cells
#include "ComplexPhiAndChiExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "CCZ4Cartoon.hpp"
#include "MovingPunctureGauge.hpp"
#include "ConstraintsCartoon.hpp"
#include "SetValue.hpp"
#include "TraceARemovalCartoon.hpp"
#include "NoetherCharge.hpp"
#include "ExpansionScalar.hpp"

// For GW extraction
#include "WeylExtraction.hpp"
#include "WeylOmScalar.hpp"

//Boson Star
#include "BosonStar.hpp"
#include "ComplexScalarField.hpp"
#include "ComplexPotential.hpp"

#include "ADMQuantities.hpp"
#include "ADMQuantitiesExtraction.hpp"
#include "GammaCartoonCalculator.hpp"

void HeadonBS2DLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance: "),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

void HeadonBS2DLevel::initialData()
{
    CH_TIME("HeadonBS2DLevel::initialData");
    if (m_verbosity)
        pout() << "HeadonBS2DLevel::initialData " << m_level << endl;

    // 1D ODE solve happens once; m_dx is patched per level before BoxLoops
    static bool s_first_call = true;
    static std::unique_ptr<BosonStar> s_boson_star;
    if (s_first_call)
    {
        pout() << "HeadonBS2DLevel::initialData on level " << m_level
               << " creating boson star profile" << endl;
        s_boson_star.reset(new BosonStar(m_p.bosonstar_params, m_p.potential_params,
                                         m_p.m_G_Newton, m_dx, m_verbosity));
        s_boson_star->compute_1d_solution(4. * m_p.L);
        s_boson_star->print_star_info();
        s_first_call = false;
    }
    s_boson_star->m_dx = m_dx; // update for this level's grid spacing

    BoxLoops::loop(make_compute_pack(SetValue(0.0), *s_boson_star), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS, disable_simd());

    BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS, disable_simd());

    fillAllGhosts();

}

// Things to do before a plot level - need to calculate the Weyl scalars
void HeadonBS2DLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    BoxLoops::loop(make_compute_pack(WeylOmScalar(m_p.extraction_params.center, m_dx),
                                    Constraints<Potential>(m_dx, potential, m_p.m_G_Newton),
                                    NoetherCharge<>(m_dx),
                                    ExpansionScalar<>(m_dx, m_p.center)), m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

void HeadonBS2DLevel::specificEvalRHS(GRLevelData &a_soln,
                                         GRLevelData &a_rhs,
                                         const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(
        make_compute_pack(TraceARemovalCartoon(), PositiveChiAndAlpha()),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate CCZ4 right hand side
    Potential potential(m_p.potential_params);
    BoxLoops::loop(
        CCZ4Cartoon<MovingPunctureGauge, FourthOrderDerivatives, Potential>(
            m_p.ccz4_params, m_dx, m_p.sigma, potential, m_p.m_G_Newton,
            m_p.formulation),
        a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

void HeadonBS2DLevel::specificUpdateODE(GRLevelData &a_soln,
                                           const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemovalCartoon(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void HeadonBS2DLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                                 const FArrayBox &current_state,
                                                 const double a_maximum)
{
    BoxLoops::loop(ComplexPhiAndChiExtractionTaggingCriterion(m_dx, m_level,
                   m_p.mass_extraction_params, m_p.regrid_threshold_phi,
                   m_p.regrid_threshold_chi), current_state, tagging_criterion);
}

void HeadonBS2DLevel::specificPostTimeStep()
{
    CH_TIME("HeadonBS2DLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main

    fillAllGhosts();

    Potential potential(m_p.potential_params);
    BoxLoops::loop(WeylOmScalar(m_p.extraction_params.center, m_dx),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);
    BoxLoops::loop(Constraints<Potential>(m_dx, potential, m_p.m_G_Newton),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);


#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
    {
        // Hack: if avg radius of found AH is negative, we reset initial guess
        double current_avg_AH_radius =  m_bh_amr.m_ah_finder.get(0)->get_ave_F();
        if (current_avg_AH_radius < 0.)
        {
            m_bh_amr.m_ah_finder.get(0)->solver.reset_initial_guess();
            pout() << "AHFinder: resetting initial guess as avg radius is "
                      "negative."
                   << endl;
        }
        else if (current_avg_AH_radius > 10.)
        {
            m_bh_amr.m_ah_finder.get(0)->solver.reset_initial_guess();
            pout() << "AHFinder: resetting initial guess as avg radius is "
                      "too large."
                   << endl;
        }
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
    }
#endif

    if (m_p.activate_extraction == 1 &&
       at_level_timestep_multiple(m_p.extraction_params.min_extraction_level()))
    {
        // Do the extraction on the min extraction level
        if (m_level == m_p.extraction_params.min_extraction_level())
        {
            if (m_verbosity)
            {
                pout() << "BinaryBSLevel::specificPostTimeStep:"
                          " Extracting gravitational waves." << endl;
            }


            // Refresh the interpolator and do the interpolation
            m_gr_amr.m_interpolator->refresh();
            WeylExtraction gw_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gw_extraction.execute_query(m_gr_amr.m_interpolator);
        }
    }

    if (at_level_timestep_multiple(0))
    {
        BoxLoops::loop(NoetherCharge<>(m_dx), m_state_new, m_state_diagnostics,
                  EXCLUDE_GHOST_CELLS);
    }

    // Robins Radius Horizon finder
    if (m_p.m_RH_activate)
    {
        // setup on the first level that fires (finest-first in GRChombo), so
        // every surface gets a t=0 entry when its own level calls update below
        if (first_step && m_bh_amr.m_rh_union.m_surfaces.empty())
            m_bh_amr.m_rh_union.setup(m_p.m_RH_num_horizons,
                                       m_p.m_RH_initial_radii,
                                       m_p.m_RH_initial_centre,
                                       m_p.m_RH_num_points,
                                       m_p.m_RH_level,
                                       m_p.m_RH_time_step_freq,
                                       m_p.m_RH_newton_crit,
                                       m_p.m_RH_chase_speeds);
        m_bh_amr.m_rh_union.update(m_time, m_level);
    }

    if (m_level == 0)
    {
        bool first_step = (m_time == 0.);
        AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
        double L2_Ham = amr_reductions.norm(c_Ham, 2, true);
        double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom2), 2, true);
        SmallDataIO constraints_file("constraint_norms",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
        constraints_file.remove_duplicate_time_data();
        if (first_step)
        {
            constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
        }
        constraints_file.write_time_data_line({L2_Ham, L2_Mom});

        int adm_min_level = 0;
        bool calculate_adm = at_level_timestep_multiple(adm_min_level);
        if (calculate_adm)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);
            double M_ADM = amr_reductions.sum(c_rho_ADM);
            SmallDataIO M_ADM_file("M_ADM", m_dt, m_time,
                                m_restart_time, SmallDataIO::APPEND,
                                first_step);
            M_ADM_file.remove_duplicate_time_data();
            if (first_step)
            {
                M_ADM_file.write_header_line({"M_ADM"});
            }
            M_ADM_file.write_time_data_line({M_ADM});
        }

        double noether_charge = amr_reductions.sum(c_N);
        SmallDataIO noether_charge_file("NoetherCharge", m_dt, m_time,
                                        m_restart_time,
                                        SmallDataIO::APPEND,
                                        first_step);
        noether_charge_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            noether_charge_file.write_header_line({"Noether Charge"});
        }
        noether_charge_file.write_time_data_line({noether_charge});

        AMRReductions<VariableType::evolution> amr_reductions_evolution(m_gr_amr);
	double min_chi = amr_reductions_evolution.min(c_chi);
        SmallDataIO min_chi_file("min_chi",
                                        m_dt, m_time, m_restart_time,
                                        SmallDataIO::APPEND, first_step);
        min_chi_file.remove_duplicate_time_data();
        if (first_step)
        {
            min_chi_file.write_header_line({"min_chi"});
        }
        min_chi_file.write_time_data_line({min_chi});

        // Compute the maximum of mod_phi and write it to a file
        double mod_phi_max = amr_reductions.max(c_mod_phi);
        SmallDataIO mod_phi_max_file("mod_phi_max", m_dt, m_time,
                                 m_restart_time,
                                 SmallDataIO::APPEND,
                                 first_step);
        mod_phi_max_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            mod_phi_max_file.write_header_line({"max mod phi"});
        }
        mod_phi_max_file.write_time_data_line({mod_phi_max});

    }

}
