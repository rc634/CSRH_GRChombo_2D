/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MOVINGBOXESREFINEMENT_HPP_
#define MOVINGBOXESREFINEMENT_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

//! This class tags cells based on three criteria - the
//! value of the second derivs, the extraction regions
//! and the puncture horizons (which must be covered to
//! a given level
class MovingBoxesRefinement
{
  protected:
    const double m_dx;
    const int m_level;
    const int m_puncture_max_level;
    const FourthOrderDerivatives m_deriv;
    const double m_puncture_radius;
    const double m_puncture_mass;
    const std::array<double, CH_SPACEDIM> m_center;
    const double m_buffer;

  public:
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };

    // The constructor
    MovingBoxesRefinement(
        const double dx, const int a_level,
        const int a_puncture_max_level,
        const std::array<double, CH_SPACEDIM> a_center,
        const double a_puncture_radius, 
        const double a_puncture_mass, 
        const double a_buffer)
        : m_dx(dx), m_level(a_level),
          m_puncture_max_level(a_puncture_max_level),
          m_deriv(dx),
          m_puncture_radius(a_puncture_radius), 
          m_puncture_mass(a_puncture_mass),
          m_center(a_center),
          m_buffer(a_buffer)
    {
    };

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // first test the gradients for regions of high curvature
        const auto d2 = m_deriv.template diff2<Vars>(current_cell);
        data_t mod_d2_chi = 0;
        FOR(idir, jdir)

         {
             mod_d2_chi += d2.chi[idir][jdir] * d2.chi[idir][jdir];
         }
        data_t criterion = m_dx * sqrt(mod_d2_chi);

        //data_t criterion;
        std::array<double, CH_SPACEDIM> puncture_centre1;

        FOR1(i) {puncture_centre1[i] =  m_center[i];}

        const Coordinates<data_t> coords(current_cell, m_dx, puncture_centre1);

		const data_t max_abs_xy =
                        simd_max(abs(coords.x), abs(coords.y));
		
		if (m_level <= m_puncture_max_level)
                    {
                        // we want the 2nd and 3rd levels above
                        // puncture_max_level to be twice the size of the next
                        // finest level
                        const double factor =
                            pow(2.0, min(m_puncture_max_level -
                                             m_level - 1, 2));

                        auto regrid = simd_compare_lt(
                            max_abs_xy,
                            factor * (m_puncture_radius / 2. +
                                      m_buffer));
			// NOTE: you can also use factor * (m_puncture_masses[ipuncture] * 2. +
                        // m_buffer) in the simd_compare_lt call; this would result in milder tagging. 
			
                        criterion = simd_conditional(regrid, 100.0, criterion);
                    }
		else
                    {
                        // remove any finer levels for BHs with
                        // puncture_max_level < max_level
                        auto dont_regrid = simd_compare_lt(
                            max_abs_xy, m_puncture_radius / 2. +
                                             m_buffer);
                        criterion =
                            simd_conditional(dont_regrid, 0.0, criterion);
                    }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* MOVINGBOXESREFINEMENT_HPP_ */
