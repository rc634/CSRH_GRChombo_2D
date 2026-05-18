/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTAR_HPP_
#define BOSONSTAR_HPP_

#include "BosonStarParams.hpp"
#include "BosonStarSolution.hpp"
#include "Cell.hpp"
#include "ComplexPotential.hpp"
//#include "ComplexScalarField.hpp"
#include "CCZ4Cartoon.hpp"
#include "Coordinates.hpp"
//#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "parstream.H" //gives pout
#include "simd.hpp"
#include <vector>

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class BosonStar
{

  public:
    //! The constructor
    BosonStar(BosonStar_params_t a_params_BosonStar,
              Potential::params_t a_params_potential, double a_G_Newton,
              double a_dx, int a_verbosity);

    //! Computes the 1d solution and stores in m_1d_sol
    void compute_1d_solution(const double max_r);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    // specific initial data computes

    // a single BH/BS
    template <class data_t> void compute_BH(Cell<data_t> current_cell) const;

    // binary spacetime BS/BH + BS/BH
    template <class data_t> void compute_BS_binary(Cell<data_t> current_cell) const;

    //! Prints ADM mass and central field amplitude to pout
    void print_star_info() const;

    // The object that stores the solution found by the 1d ODE integrator
    BosonStarSolution m_1d_sol;
    double m_dx;

  protected:
    double m_G_Newton;
    BosonStar_params_t m_params_BosonStar;
    Potential::params_t m_params_potential; //!< The potential params
    int m_verbosity;
};

#include "BosonStar.impl.hpp"

#endif /* BOSONSTAR_HPP_ */
