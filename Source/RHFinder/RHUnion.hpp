#ifndef RHUNION_HPP_
#define RHUNION_HPP_

#include "AMRInterpolator.hpp"
#include "InterpolationQuery.hpp"
#include "Lagrange.hpp"
#include "RHSurf.hpp"
#include "UserVariables.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

// Owns all tracked horizons.
class RHUnion
{
  public:
    enum : int { NG = 2 }; // ghost cells on each side of each surface array

    int    m_n_chase           = 100;     // max interpolation calls per update
    int    m_num_stale_repeats = 5;       // stale repeats when err > thresh_high
    double m_courant           = 0.25;     // normal chase step size
    double m_thresh_high       = 0.03;    // above: slow chase + stale repeats
    double m_thresh_close     = 0.0005;   // below: switch from fast chase to Newton
    double m_thresh_super_low  = 0.000001; // below: surface converged, stop
    double m_newton_delta_f    = 1e-4;    // finite difference step for Jacobian assembly
    std::vector<RHSurf> m_surfaces;
    std::vector<std::ofstream> m_outfiles;

    RHUnion() {}

    void setup(int a_num,
               const std::vector<double> &a_radii,
               const std::vector<double> &a_centres_x,
               const std::vector<int>    &a_num_points,
               const std::vector<int>    &a_levels,
               const std::vector<int>    &a_time_step_freqs,
               const std::vector<double> &a_newton_crits)
    {
        m_surfaces.clear();
        m_outfiles.clear();
        m_outfiles.resize(a_num);
        for (int i = 0; i < a_num; ++i)
        {
            m_surfaces.emplace_back(a_num_points[i], NG, std::vector<double>{a_centres_x[i], 0.0}, i);
            m_surfaces.back().m_level          = a_levels[i];
            m_surfaces.back().m_time_step_freq = a_time_step_freqs[i];
            m_surfaces.back().m_newton_crit    = a_newton_crits[i];
            std::fill(m_surfaces.back().m_f.begin(),
                      m_surfaces.back().m_f.end(), a_radii[i]);
            if (procID() == 0)
            {
                m_outfiles[i].open("rh_surf_" + std::to_string(i) + ".dat",
                                   std::ios::out | std::ios::trunc);
                m_outfiles[i] << std::setw(16) << "#t"
                              << std::setw(6)  << "surf"
                              << std::setw(16) << "x"
                              << std::setw(16) << "<r>"
                              << std::setw(16) << "Area"
                              << std::setw(16) << "M_irr"
                              << std::setw(16) << "P_x"
                              << std::setw(16) << "<Theta+>"
                              << std::setw(16) << "<Theta->"
                              << std::setw(16) << "err"
                              << std::setw(8)  << "mode"
                              << std::endl;
            }
        }

        if (procID() == 0)
        {
            pout() << "\nRHUnion::Construction  Initial surfaces are:" << std::endl;
            for (const auto &surf : m_surfaces)
                surf.hello();
        }
    }

    void set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        m_interpolator = a_interpolator;
    }

    // Interpolate all needed fields at interior surface points in one combined
    // query, then store results in each surface's arrays and fill ghost cells.
    void interpolate_fields()
    {
        int total_pts = 0;
        for (const auto &surf : m_surfaces) {
            total_pts += surf.m_n;
        }

        std::vector<double> qx(total_pts), qy(total_pts);
        int offset = 0;
        for (const auto &surf : m_surfaces)
        {
            for (int i = 0; i < surf.m_n; ++i)
            {
                int ii = i + surf.m_NG; // non-ghost array index
                qx[offset + i] = surf.m_centre[0] + surf.m_f[ii] * std::cos(surf.m_theta[ii]);
                qy[offset + i] = surf.m_centre[1] + surf.m_f[ii] * std::sin(surf.m_theta[ii]);
            }
            offset += surf.m_n;
        }

        // output buffers — values
        std::vector<double> chi(total_pts), K(total_pts);
        std::vector<double> A11(total_pts), A12(total_pts), A22(total_pts), Aww(total_pts);
        std::vector<double> h11(total_pts), h12(total_pts), h22(total_pts), hww(total_pts);
        // output buffers — derivatives of chi and conformal metric
        std::vector<double> dx_chi(total_pts), dy_chi(total_pts);
        std::vector<double> dx_h11(total_pts), dy_h11(total_pts);
        std::vector<double> dx_h12(total_pts), dy_h12(total_pts);
        std::vector<double> dx_h22(total_pts), dy_h22(total_pts);
        std::vector<double> dx_hww(total_pts), dy_hww(total_pts);

        InterpolationQuery query(total_pts);
        query.setCoords(0, qx.data()).setCoords(1, qy.data());

        query.addComp(c_chi, chi.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_K,   K.data(),   Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_A11, A11.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_A12, A12.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_A22, A22.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_Aww, Aww.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_h11, h11.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_h12, h12.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_h22, h22.data(), Derivative::LOCAL, VariableType::evolution);
        query.addComp(c_hww, hww.data(), Derivative::LOCAL, VariableType::evolution);

        query.addComp(c_chi, dx_chi.data(), Derivative::dx, VariableType::evolution);
        query.addComp(c_chi, dy_chi.data(), Derivative::dy, VariableType::evolution);
        query.addComp(c_h11, dx_h11.data(), Derivative::dx, VariableType::evolution);
        query.addComp(c_h11, dy_h11.data(), Derivative::dy, VariableType::evolution);
        query.addComp(c_h12, dx_h12.data(), Derivative::dx, VariableType::evolution);
        query.addComp(c_h12, dy_h12.data(), Derivative::dy, VariableType::evolution);
        query.addComp(c_h22, dx_h22.data(), Derivative::dx, VariableType::evolution);
        query.addComp(c_h22, dy_h22.data(), Derivative::dy, VariableType::evolution);
        query.addComp(c_hww, dx_hww.data(), Derivative::dx, VariableType::evolution);
        query.addComp(c_hww, dy_hww.data(), Derivative::dy, VariableType::evolution);

        m_interpolator->interp(query);

        offset = 0;
        for (auto &surf : m_surfaces)
        {
            for (int i = 0; i < surf.m_n; ++i)
            {
                int ii = surf.m_NG + i; // non-ghost array index
                surf.m_chi[ii]    = chi[offset + i];
                surf.m_K[ii]      = K[offset + i];
                surf.m_A11[ii]    = A11[offset + i];
                surf.m_A12[ii]    = A12[offset + i];
                surf.m_A22[ii]    = A22[offset + i];
                surf.m_Aww[ii]    = Aww[offset + i];
                surf.m_h11[ii]    = h11[offset + i];
                surf.m_h12[ii]    = h12[offset + i];
                surf.m_h22[ii]    = h22[offset + i];
                surf.m_hww[ii]    = hww[offset + i];
                surf.m_dx_chi[ii] = dx_chi[offset + i];
                surf.m_dy_chi[ii] = dy_chi[offset + i];
                surf.m_dx_h11[ii] = dx_h11[offset + i];
                surf.m_dy_h11[ii] = dy_h11[offset + i];
                surf.m_dx_h12[ii] = dx_h12[offset + i];
                surf.m_dy_h12[ii] = dy_h12[offset + i];
                surf.m_dx_h22[ii] = dx_h22[offset + i];
                surf.m_dy_h22[ii] = dy_h22[offset + i];
                surf.m_dx_hww[ii] = dx_hww[offset + i];
                surf.m_dy_hww[ii] = dy_hww[offset + i];
            }
            surf.fill_all_ghosts();
            offset += surf.m_n;
        }
    }

    void update(double a_time, int a_level)
    {
        // skip entirely if no surface runs on this AMR level
        bool any = false;
        for (const auto &s : m_surfaces)
            if (s.m_level == a_level) { any = true; break; }
        if (!any) return;

        m_interpolator->refresh();
        interpolate_fields();

        for (auto &surf : m_surfaces)
            surf.re_centre();

        const int n = (int)m_surfaces.size();

        // per-surface error and state (all surfaces, so diagnostics stay current)
        std::vector<double> errs(n);
        for (int k = 0; k < n; ++k)
        {
            errs[k] = m_surfaces[k].expansion_error();
            auto &surf = m_surfaces[k];
            if (errs[k] <= m_thresh_super_low)
                surf.m_state = RHSurf::SolverState::FOUND;
            else if (errs[k] <= m_thresh_close)
                surf.m_state = RHSurf::SolverState::CLOSE;
            else if (errs[k] <= m_thresh_high)
                surf.m_state = RHSurf::SolverState::MEDIUM;
            else
                surf.m_state = RHSurf::SolverState::FAR;
        }

        // chase + Newton, per surface
        for (int k = 0; k < n; ++k)
        {
            auto &surf = m_surfaces[k];
            if (surf.m_level != a_level) continue;
            if (surf.m_dead) continue;
            if (errs[k] <= m_thresh_super_low) continue;

            try
            {
                const double err_pre = errs[k];

                if (err_pre > m_thresh_high)
                {
                    for (int i = 0; i < surf.m_time_step_freq; ++i)
                    {
                        surf.chase_step(m_courant);
                        for (int j = 0; j < m_num_stale_repeats; ++j)
                            surf.chase_step(m_courant);
                        interpolate_fields();
                    }
                }
                else
                {
                    for (int i = 0; i < surf.m_time_step_freq; ++i)
                    {
                        surf.chase_step(m_courant);
                        interpolate_fields();
                    }
                }

                if (surf.m_newton_crit > 0.0)
                {
                    const double err_post = surf.expansion_error();
                    if (err_post < surf.m_newton_crit)
                    {
                        const std::vector<double> f_before_newton = surf.m_f;

                        for (int ns = 0; ns < 10; ++ns)
                        {
                            interpolate_fields();
                            surf.banded_newton_step(m_newton_delta_f, m_courant);
                        }
                        interpolate_fields();

                        // Revert if Newton made things worse
                        const double err_after_newton = surf.expansion_error();
                        if (err_after_newton > err_post)
                        {
                            surf.m_f = f_before_newton;
                            surf.fill_all_ghosts();
                            interpolate_fields();
                            pout() << "RHFinder: Newton worsened surface " << k
                                   << " (" << err_post << " -> " << err_after_newton
                                   << "), reverting.\n";
                        }
                    }
                }
            }
            catch (const std::exception &e)
            {
                surf.m_dead = true;
                pout() << "\nRHFinder: surface " << k << " marked dead: "
                       << e.what() << "\n";
            }
        }


        // re-evaluate states so fancy_hello and rhout reflect post-solve error
        for (int k = 0; k < n; ++k)
        {
            auto &surf = m_surfaces[k];
            if (surf.m_dead) continue;
            const double err = surf.expansion_error();
            if (err <= m_thresh_super_low)   surf.m_state = RHSurf::SolverState::FOUND;
            else if (err <= m_thresh_close) surf.m_state = RHSurf::SolverState::CLOSE;
            else if (err <= m_thresh_high)   surf.m_state = RHSurf::SolverState::MEDIUM;
            else                             surf.m_state = RHSurf::SolverState::FAR;
        }

        // file outputs, rh_out ...
        rhout(a_time);

        if (procID() == 0)
        {
            pout() << "\nRHFinder::update  t = " << std::fixed << std::setprecision(4) << a_time;
            for (const auto &surf : m_surfaces)
                surf.fancy_hello();
        }
    }

    void rhout(double a_time)
    {
        if (procID() != 0)
        {
            return;
        }
        for (int i = 0; i < (int)m_surfaces.size(); ++i)
        {
            m_surfaces[i].print_diagnostics(m_outfiles[i], a_time);
        }
    }

    void hello() const
    {
        for (const auto &surf : m_surfaces)
        {
            surf.hello();
        }
    }

  private:
    AMRInterpolator<Lagrange<4>> *m_interpolator = nullptr;
};

#endif /* RHUNION_HPP_ */
