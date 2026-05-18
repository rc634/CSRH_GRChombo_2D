/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
#include "ComplexPotential.hpp"

// Problem specific(from robin) includes:
#include "RotatingBosonStarParams.hpp"
//#include "ComplexPotential.hpp"

#ifdef USE_AHFINDER
#include "AHInitialGuess.hpp"
#endif

class SimulationParameters : public SimulationParametersBase
{
public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        readParams(pp);
        check_params();
    }

    void readParams(GRParmParse &pp)
    {
        pp.load("G_Newton", m_G_Newton, 1.0);

        pp.load("star_centre", rotating_bosonstar_params.star_centre,
                {0, 0});
        pp.load("initial_data_path", rotating_bosonstar_params.base_path, std::string(""));
        pp.load("BS_frequency", rotating_bosonstar_params.BS_frequency, 0.15910835770266477);

        //Tagging
        pp.load("star_radius", puncture_radius, 4.);
        pp.load("star_mass", puncture_mass, 4.);
        pp.load("tag_buffer", tag_buffer, 0.5);
	pp.load("tag_puncture_max_level", tag_puncture_max_level, max_level);

        // Potential params
        pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
        pp.load("phi4_coeff", potential_params.phi4_coeff, 0.0);
        pp.load("solitonic", potential_params.solitonic, false);
        pp.load("sigma_soliton", potential_params.sigma_soliton, 0.02);

        // BubbleTaggingCriterion regridding
        pp.load("threshold_phi", regrid_threshold_phi, 1.);
        pp.load("threshold_chi", regrid_threshold_chi, 1.); 

        // Do we want Weyl extraction, puncture tracking and constraint norm
        // calculation?
        pp.load("activate_extraction", activate_extraction, false);

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("num_mass_extraction_radii",
                mass_extraction_params.num_extraction_radii, 1);
        pp.load("mass_extraction_levels",
                mass_extraction_params.extraction_levels,
                mass_extraction_params.num_extraction_radii, 0);
        pp.load("mass_extraction_radii",
                mass_extraction_params.extraction_radii,
                mass_extraction_params.num_extraction_radii, 0.1);
        pp.load("num_points_phi_mass", mass_extraction_params.num_points_phi,
                2);
        pp.load("num_points_theta_mass",
                mass_extraction_params.num_points_theta, 4);
        pp.load("mass_extraction_center",
                mass_extraction_params.extraction_center,
                {0.0, 0.0});

#ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess,
                0.5 * puncture_mass);
#endif
    }

    void check_params()
    {}

    double m_G_Newton;
    bool activate_extraction;

    RotatingBosonStar_params_t rotating_bosonstar_params;
    Potential::params_t potential_params;

    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;
    bool identical;
    
    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi;
    Real tag_buffer;
    Real puncture_mass, puncture_radius;
    int tag_puncture_max_level;
    
#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif

#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif
};
#endif /* SIMULATIONPARAMETERS_HPP_ */
