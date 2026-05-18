/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef THINSHELL2DLEVEL_HPP_
#define THINSHELL2DLEVEL_HPP_

#include "STAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
#include "ComplexPotential.hpp"
#include "ComplexScalarField.hpp"
#include "STAMR.hpp"

class ThinShell2DLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<ThinShell2DLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    STAMR &m_st_amr = dynamic_cast<STAMR &>(m_gr_amr);

     // Typedef for scalar field
    typedef ComplexScalarField<Potential> ComplexScalarFieldWithPotential;

    /// Things to do at every full timestep
    ///(might include several substeps, e.g. in RK4)
    virtual void specificAdvance() override;

    /// Initial data calculation
    virtual void initialData() override;

    /// Any actions that should happen just before plot files output
    virtual void prePlotLevel() override;

    /// Calculation of the right hand side for the time stepping
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state,
                                         const double a_maximum) override;

    virtual void specificPostTimeStep() override;
};

#endif /* THINSHELL2DLEVEL_HPP_ */
