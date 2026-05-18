/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EFFECTIVEPOTENTIALEXTRACTION_HPP_
#define EFFECTIVEPOTENTIALEXTRACTION_HPP_

#include "SphericalExtraction.hpp"
#include "simd.hpp"

class EffectivePotentialExtraction : public SphericalExtraction
{
  public:
    //! The constructor
    EffectivePotentialExtraction(const spherical_extraction_params_t &a_params, double a_dt,
                   double a_time, bool a_first_step,
                   double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        add_var(c_proper_dist, VariableType::diagnostic);
        add_var(c_lapse_sq, VariableType::diagnostic);
        add_var(c_beta_sq, VariableType::diagnostic);
    }

    void run(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        std::vector<double> integrals1;
        std::vector<double> integrals2;
        std::vector<double> integrals3;

        auto integrand1 = [](std::vector<double> diagnostic_vars, double r,
                            double theta,
                            double phi) { return diagnostic_vars[0]; };

        auto integrand2 = [](std::vector<double> diagnostic_vars, double r,
                            double theta,
                            double phi) { return diagnostic_vars[2]; };
        
        auto integrand3 = [](std::vector<double> diagnostic_vars, double r,
                            double theta,
                            double phi) { return diagnostic_vars[1]; };

        extract(a_interpolator);
        add_integrand(integrand1, integrals1);
        integrate();

        pout() << "G_phi_phi " << integrals1[1] << endl;

        add_integrand(integrand2, integrals2);
        integrate();

        pout() << "Beta_sq " << integrals2[1] << endl;

        add_integrand(integrand3, integrals3);
        integrate();

        std::vector<double> lapse_sq = integrals3;

        std::transform(integrals1.begin(), integrals1.end(), integrals1.begin(),
               [](double value) { return (value * value / (4 * M_PI * M_PI)); });
        
        std::transform(integrals2.begin(), integrals2.end(), integrals2.begin(),
               [](double value) { return (value / (2 * M_PI)); });

        std::transform(integrals3.begin(), integrals3.end(), integrals3.begin(),
               [](double value) { return (value / (2 * M_PI)); });

        std::vector<double> eff_potential;
        eff_potential.resize(integrals1.size());
        for (int i = 0; i < integrals1.size(); ++i)
        {
            eff_potential[i] = (sqrt(integrals1[i]) * sqrt(integrals2[i]) - sqrt(integrals1[i] * integrals2[i] - (-integrals3[i] + integrals2[i]) * integrals1[i])) / (integrals1[i]);
        }

        std::string integrals_filename = "eff_potential";
        std::vector<std::vector<double>> integrals_for_writing = {eff_potential};
        std::vector<std::string> labels = {"V_eff"};
        write_integrals(integrals_filename, integrals_for_writing, labels);
    }
};

#endif /*EFFECTIVEPOTENTIALEXTRACTION_HPP_ */