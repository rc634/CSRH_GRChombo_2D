#pragma once

#include <vector>   // add this
#include <cmath>




#ifndef BSPIXEL_HPP_
#define BSPIXEL_HPP_
#include "BosonStarSolution.hpp"
#include "BosonStarParams.hpp"

struct BSPixel
{
    // Scalar field
    double phi    = 0.;
    double phi_Im = 0.;
    double Pi     = 0.;
    double Pi_Im  = 0.;

    // Spatial metric (lower indices)
    double chi = 1.;
    double gamma[3][3] = {
        {1., 0., 0.},
        {0., 1., 0.},
        {0., 0., 1.}
    };

    // Extrinsic curvature (lower indices, symmetric)
    double trK = 0.;
    double K[3][3] = {
        {0., 0., 0.},
        {0., 0., 0.},
        {0., 0., 0.}
    };

    // Gauge
    double lapse   = 1.;
    double shift_x = 0.;

    // other
    double w      = 0.;
    double direction = 0.; // -1 = left
    double m = 0.;
    bool is_BH = false; //


    void init(double a_dir,
              Coordinates<double> coords,
              const BosonStarSolution& sol,
              const BosonStar_params_t& BSpp) {
        // - * - * - * - * - * - * - * - * //
        // load everything about a single BS that the pixel needs to see


        // params
        w = sol.get_w();
        if (is_BH && BSpp.BH_override)
        {
            m = BSpp.BH_mass;
        }
        else
        {
            m = sol.get_ADMmass();
        }


        // Boosts
        direction = a_dir;
        const double rapidity = direction * BSpp.binary_rapidity;
        const double c_ = cosh(rapidity);
        const double s_ = sinh(rapidity);
        const double v_ = tanh(rapidity);


        // Coordinates
        // * - * - * -
        // lab coords -tilde-
        const double x_tilde = coords.x + 0.5 * direction * BSpp.binary_separation;
        const double t_tilde = 0.;
        // rest frame coords
        const double x = c_ * x_tilde + s_ * t_tilde;
        const double y = coords.y;
        const double z = 0.;
        const double t = c_ * t_tilde - s_ * x_tilde;
        const double r = sqrt(x*x + y*y);
        const double safe_r = sqrt(r*r + 10e-16);


        // fields
        double p = sol.get_p_interp(r);
        double dp = sol.get_dp_interp(r);
        double omega = sol.get_lapse_interp(r);
        double domega = sol.get_dlapse_interp(r);
        double psi = sol.get_psi_interp(r);
        double dpsi = sol.get_dpsi_interp(r);
        double pc_os = psi*psi*c_*c_ - omega*omega*s_*s_;
        if (is_BH) {
            p = 0.;
            dp = 0.;
            psi = pow(1.+m/(2.*safe_r),2);
            dpsi = -(1.+m/(2.*safe_r))*m/(safe_r*safe_r);
            omega = (2.*r-m)/(2.*r+m);
            domega = 4.*m / pow(2.*r+m, 2);
            pc_os = psi*psi*c_*c_ - omega*omega*s_*s_;
        }


        // analytic gauge vars
        lapse = omega * psi / sqrt(pc_os);
        shift_x = s_ * c_ * (psi*psi - omega*omega) / pc_os;


        // scalar feild output vars
        const double phase = BSpp.BS_phase * M_PI + w * t;
        phi    = p * cos(phase);
        phi_Im = p * sin(phase);
        Pi    = -(1./lapse) * (
                    (x_tilde/r) * (s_ - shift_x*c_) * dp * cos(phase)
                               - w * (c_ - shift_x*s_) * p * sin(phase) );
        Pi_Im = -(1./lapse) * (
                    (x_tilde/r) * (s_ - shift_x*c_) * dp * sin(phase)
                               + w * (c_ - shift_x*s_) * p * cos(phase) );


        // spatial metric
        gamma[0][0] = pc_os;
        gamma[1][1] = psi * psi;
        gamma[2][2] = psi * psi;
        chi = pow(gamma[0][0]*gamma[1][1]*gamma[2][2],-1./3.);


        // extrinsic curvature
        const double dpsi_o_psi     = dpsi/psi;
        const double domega_o_omega = domega/omega;
        K[2][2] = lapse * s_ * (x_tilde/r) * dpsi_o_psi;
        K[1][1] = K[2][2];
        K[0][1] = lapse * c_ * s_ * (y/r) * (domega_o_omega - dpsi_o_psi);
        K[0][2] = lapse * c_ * s_ * (z/r) * (domega_o_omega - dpsi_o_psi);
        K[1][2] = 0.;
        K[1][0] = K[0][1];
        K[2][0] = K[0][2];
        K[2][1] = K[1][2];
        K[0][0] = lapse * (x_tilde/r) * s_ * c_* c_
                   * ( 2.*domega_o_omega - dpsi_o_psi
                      - v_*v_ * omega * domega / (psi*psi));

        // TRACE K calculated with DIAGONAL METRIC
        for (size_t i = 0; i < 3; i++) {
          trK += K[i][i] / gamma[i][i];
        }


        // end * end * end * end * end * ------ //
    }
};


#endif /* BSPIXEL_HPP_ */
