

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
                                = CCZ4CartoonVars::VarsWithGauge<data_t>;;

public:

    NoetherCharge(const double a_dx) : m_deriv(a_dx) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // load vars locally
        const auto vars = current_cell.template load_vars<Vars>();
        // const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
        // const auto matter_vars = current_cell.template load_vars<MatterVars>();
        // const auto advec_csf =
        //     this->m_deriv.template advection<Vars>(current_cell, vars.shift);

        // Coordinates<data_t> coords(current_cell, this->m_deriv.m_dx);
        //
        // using namespace TensorAlgebra;

        // calculate Noether charge
        data_t N = pow(vars.chi, -1.5) * (vars.phi_Im
            * vars.Pi - vars.phi * vars.Pi_Im);

        //Calculate modulus of the scalar field
        data_t mod_phi = sqrt(vars.phi * vars.phi
                            + vars.phi_Im * vars.phi_Im);

        current_cell.store_vars(N, c_N);
        current_cell.store_vars(mod_phi, c_mod_phi);
    }
};

#endif /* NOETHERCHARGE_HPP_ */
