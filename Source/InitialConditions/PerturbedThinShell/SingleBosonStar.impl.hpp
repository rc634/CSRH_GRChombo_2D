/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SINGLEBOSONSTAR_HPP_)
#error "This file should only be included through SingleBosonStar.hpp"
#endif

#ifndef SINGLEBOSONSTAR_IMPL_HPP_
#define SINGLEBOSONSTAR_IMPL_HPP_

#include <iostream>
#include <vector>
#include <fstream> // Needed for std::ifstream

SingleBosonStar::SingleBosonStar() {}

// Read below 
void SingleBosonStar::main(std::string m_base_path)
{
    std::vector<double> temp_data_x1;  // array for theta & array for radius 
    double value_x1;
    // std::string base_path = "/Users/macbookpro/Desktop/PHD/Projects/InstabilityRing/";
    std::string filename_radius = "radius.dat";
    std::string filename_amp = "phi_Re.dat";
    std::string filename_conformal_factor = "chi.dat";
    std::string filename_lapse = "lapse.dat";
    std::string filename_pi = "Pi_Im.dat";

    std::string full_path_radius = m_base_path + filename_radius;
    std::string full_path_amp = m_base_path + filename_amp;
    std::string full_path_conformal_factor = m_base_path + filename_conformal_factor;
    std::string full_path_lapse = m_base_path + filename_lapse;
    std::string full_path_pi = m_base_path + filename_pi;

    ////////////
    // Radius //
    ////////////

    std::ifstream file_x1(full_path_radius);
    if (file_x1.is_open()) {
        // Read each value from the file
        while (file_x1 >> value_x1) {
            temp_data_x1.push_back(value_x1);  // Append to temporary vector
        }
        file_x1.close();
    } else {
        std::cerr << "Unable to open file radius.dat" << std::endl;
    }
    m = temp_data_x1.size();
    x1.resize(m);
    for (int i = 0; i < m; ++i) {
        x1[i] = temp_data_x1[i];  // Copy from std::vector to VecDoub
    }

    // Resize all matrix variables
    amp.resize(m);
    chi.resize(m);     
    lapse.resize(m);
    pi_im.resize(m);    

    ///////////////
    // Amplitude //
    ///////////////

    std::ifstream file_A(full_path_amp);
    if (file_A.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < m; ++i) {
                file_A >> amp[i];
            }
        file_A.close();
    } else {
            std::cerr << "Unable to open file phi_Re.dat" << std::endl;
        }

    ////////////
    // lapse //
    ///////////

    std::ifstream file_l(full_path_lapse);
    if (file_l.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < m; ++i) {
                file_l >> lapse[i];
        }
        file_l.close();
        } else {
            std::cerr << "Unable to open file lapse.dat" << std::endl;
        }
    
    //////////////////////
    // conformal factor //
    /////////////////////

    std::ifstream file_chi(full_path_conformal_factor);
    if (file_chi.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < m; ++i) {
                file_chi >> chi[i];
        }
        file_chi.close();
        } else {
            std::cerr << "Unable to open file chi.dat" << std::endl;
        }

    //////////////////////
    //   imaginary pi   //
    /////////////////////

    std::ifstream file_pi(full_path_pi);
    if (file_pi.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < m; ++i) {
                file_pi >> pi_im[i];
        }
        file_pi.close();
        } else {
            std::cerr << "Unable to open file Pi_Im.dat" << std::endl;
        }

}

double SingleBosonStar::get_amp_interp(double r, int m) const
{
    double interpolated_value; 
    if (x1.size()==0) {
    MayDay::Error("Error: x1 is empty. Make sure to initialize it before calling get_amp_interp().");
    }

    Spline_interp cubic_inerpolation(x1,amp); 
    interpolated_value = cubic_inerpolation.interp(r);
    
    return interpolated_value;
}

double SingleBosonStar::get_pi_interp(double r, int m) const
{
    double interpolated_value; 
    if (x1.size()==0) {
    MayDay::Error("Error: x1 is empty. Make sure to initialize it before calling get_pi_interp().");
    }

    Spline_interp cubic_inerpolation(x1,pi_im); 
    interpolated_value = cubic_inerpolation.interp(r);
    
    return interpolated_value;
}

double SingleBosonStar::get_conformal_factor_interp(const double r, int m) const
{
    double interpolated_value; 

    Spline_interp cubic_inerpolation(x1,chi); 
    interpolated_value = cubic_inerpolation.interp(r);
    
    return interpolated_value;
}

double SingleBosonStar::get_lapse_interp(const double r, int m) const
{
    double interpolated_value; 

    Spline_interp cubic_inerpolation(x1,lapse); 
    interpolated_value = cubic_inerpolation.interp(r);
    
    return interpolated_value;
}

#endif /* SINGLEBOSONSTAR_IMPL_HPP_ */
