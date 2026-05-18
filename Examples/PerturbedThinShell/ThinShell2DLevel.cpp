/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */
#include <fstream>
#include <mutex>
#include <memory>

#include "ThinShell2DLevel.hpp"

// General headers
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"

// For tag cells
#include "ComplexPhiAndChiExtractionTaggingCriterion.hpp"
#include "MovingBoxesRefinement.hpp"

// Problem specific includes
#include "CCZ4Cartoon.hpp"
#include "MovingPunctureGauge.hpp"
#include "ConstraintsCartoon.hpp"
#include "SetValue.hpp"
#include "TraceARemovalCartoon.hpp"
#include "NoetherCharge.hpp"
#include "DiagnosticAij.hpp" 

// For GW extraction
#include "WeylExtraction.hpp"
#include "WeylOmScalar.hpp"

//Boson Star
#include "ThinShellInitialData.hpp"
#include "ComplexScalarField.hpp"
#include "ComplexPotential.hpp"

// ADM vars 
#include "ADMQuantities.hpp"
#include "ADMQuantitiesExtraction.hpp"
#include "GammaCartoonCalculator.hpp"

void ThinShell2DLevel::specificAdvance()
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

void ThinShell2DLevel::initialData()
{
    CH_TIME("ThinShell2DLevel::initialData");
    if (m_verbosity)
        pout() << "ThinShell2DLevel::initialData " << m_level << endl;

    // When changing class here, don't forget to change potential if necessary
    ThinShellInitialData thin_shell_initial_data(m_p.thin_shell_params, m_dx);

    if (m_verbosity)
        pout() << "BosonStarLevel::initialData & compute_1d_solution " << m_level << endl;

    thin_shell_initial_data.compute_1d_solution();

    BoxLoops::loop(make_compute_pack(SetValue(0.0), thin_shell_initial_data), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS, disable_simd());

    BoxLoops::loop(GammaCartoonCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS, disable_simd());

    fillAllGhosts(); 
}

// Things to do before a plot level - need to calculate the Weyl scalars
void ThinShell2DLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    // ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(make_compute_pack(WeylOmScalar(m_p.extraction_params.center, m_dx), 
                                    Constraints<Potential>(m_dx, potential, m_p.m_G_Newton), 
                                    NoetherCharge()), m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

void ThinShell2DLevel::specificEvalRHS(GRLevelData &a_soln,
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

void ThinShell2DLevel::specificUpdateODE(GRLevelData &a_soln,
                                           const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemovalCartoon(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void ThinShell2DLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                                 const FArrayBox &current_state,
                                                 const double a_maximum)
{
    // BoxLoops::loop(ComplexPhiAndChiExtractionTaggingCriterion(m_dx, m_level,
    //                m_p.mass_extraction_params, m_p.regrid_threshold_phi,
    //                m_p.regrid_threshold_chi, m_p.activate_extraction), current_state, tagging_criterion);
    const std::vector<double> star_coords =
            m_st_amr.m_star_tracker.get_puncture_coords();

    BoxLoops::loop(MovingBoxesRefinement(
                           m_dx, m_level, m_p.tag_puncture_max_level,
                           star_coords, m_p.do_star_track, m_p.puncture_radius, m_p.puncture_mass, m_p.tag_buffer),
                      current_state, tagging_criterion);	
}

void ThinShell2DLevel::specificPostTimeStep()
{
    CH_TIME("ThinShell2DLevel::specificPostTimeStep");

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

    if (m_p.do_star_track && m_level == m_p.star_track_level)
    {
        pout() << "Running a star tracker now" << endl;
        // if at restart time read data from dat file,
        // will default to param file if restart time is 0
        if (fabs(m_time - m_restart_time) < m_dt * 1.1)
        {
            m_st_amr.m_star_tracker.read_old_centre_from_dat(
                "StarCentres", m_dt, m_time, m_restart_time, first_step);
        }
        m_st_amr.m_star_tracker.update_star_centres(m_dt);
        m_st_amr.m_star_tracker.write_to_dat("StarCentres", m_dt, m_time,
                                             m_restart_time, first_step);
    }

#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
    {
        m_st_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
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
        BoxLoops::loop(NoetherCharge(), m_state_new, m_state_diagnostics,
                  EXCLUDE_GHOST_CELLS);
        BoxLoops::loop(DiagnosticAij<MovingPunctureGauge, FourthOrderDerivatives, Potential>(
            m_p.ccz4_params, m_dx, m_p.sigma, potential, m_p.m_G_Newton,
            m_p.formulation), m_state_new, m_state_diagnostics,
                  EXCLUDE_GHOST_CELLS);
    }

    if (m_level == 0)
    {
        bool first_step = (m_time == 0.);
        AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
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
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
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

        double A11_max = amr_reductions.max(c_dA11);
        double A22_max = amr_reductions.max(c_dA22);
        double A12_max = amr_reductions.max(c_dA12);
        double Aww_max = amr_reductions.max(c_dAww);
        double K_max = amr_reductions.max(c_dK);

        double A11_min = amr_reductions.min(c_dA11);
        double A22_min = amr_reductions.min(c_dA22);
        double A12_min = amr_reductions.min(c_dA12);
        double Aww_min = amr_reductions.min(c_dAww);
        double K_min = amr_reductions.min(c_dK);

        SmallDataIO spherical_diagnostics("SphericalDiag", m_dt, m_time,
                                        m_restart_time,
                                        SmallDataIO::APPEND,
                                        first_step);
        spherical_diagnostics.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            spherical_diagnostics.write_header_line({"A11max", "A22max", "A12max", "Awwmax", "Kmax", "A11min", "A22min", "A12min", "Awwmin", "Kmin"});
        }
        spherical_diagnostics.write_time_data_line({A11_max, A22_max, A12_max, Aww_max, K_max, A11_min, A22_min, A12_min, Aww_min, K_min});

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
