/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EFFECTIVEPOTENTIAL_HPP_
#define EFFECTIVEPOTENTIAL_HPP_

#include "CCZ4CartoonVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"
#include <array>

#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

class EffectivePotential
{
  public:
    // Use the variable definitions containing the needed quantities
    template <class data_t> using Vars = CCZ4CartoonVars::VarsWithGauge<data_t>;
    template <class data_t>
    using Diff2Vars = CCZ4CartoonVars::Diff2VarsNoGauge<data_t>;

    EffectivePotential(const std::array<double, CH_SPACEDIM> a_center,
                 const double a_dx)
        : m_center(a_center), m_deriv(a_dx)
    {
    }

    //! The compute member which calculates the wave quantities at each point on
    //! the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center
    const FourthOrderDerivatives m_deriv; //!< for calculating derivs of vars
};

#include "EffectivePotential.impl.hpp"

#endif /* EFFECTIVEPOTENTIAL_HPP_ */