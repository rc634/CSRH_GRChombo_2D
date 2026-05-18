/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ROTATINGBOSONSTAR_HPP_
#define ROTATINGBOSONSTAR_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
// #include "MatterCCZ4.hpp"
#include "ComplexScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "RotatingBosonStarSolution.hpp"
#include "RotatingBosonStarParams.hpp"
#include "ComplexPotential.hpp"
#include <vector>
#include "parstream.H" //gives pout

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class RotatingBosonStar
{

public:
    //! The constructor
    RotatingBosonStar(RotatingBosonStar_params_t a_params_RotatingBosonStar, double a_dx);

    //! Computes the 1d spinning solution 
    void compute_1d_rotating_solution();

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t>
    void compute(Cell<data_t> current_cell) const;

    RotatingBosonStarSolution rotating_BS_sol;

protected:
    RotatingBosonStar_params_t m_params_RotatingBosonStar;
    double m_dx;
};

#include "RotatingBosonStar.impl.hpp"

#endif /* ROTATINGBOSONSTAR_HPP_ */
