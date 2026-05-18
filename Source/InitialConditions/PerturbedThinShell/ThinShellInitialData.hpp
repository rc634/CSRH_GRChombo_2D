/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef THINSHELLINITIALDATA_HPP_
#define THINSHELLINITIALDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "ComplexScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "SingleBosonStar.hpp"
#include "ThinShellParams.hpp"
#include "ComplexPotential.hpp"
#include <vector>
#include "parstream.H" //gives pout

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class ThinShellInitialData
{

public:
    //! The constructor
    ThinShellInitialData(ThinShell_params_t a_params_ThinShell, double a_dx);

    //! Computes the 1d spinning solution 
    void compute_1d_solution();

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t>
    void compute(Cell<data_t> current_cell) const;

    SingleBosonStar thin_shell_sol;

protected:
    ThinShell_params_t m_params_ThinShell;
    double m_dx;
};

#include "ThinShellInitialData.impl.hpp"

#endif /* THINSHELLINITIALDATA_HPP_ */