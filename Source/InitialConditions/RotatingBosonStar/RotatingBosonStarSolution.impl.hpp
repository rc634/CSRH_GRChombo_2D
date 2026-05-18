/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(ROTATINGBOSONSTARSOLUTION_HPP_)
#error "This file should only be included through RotatingBosonStarSolution.hpp"
#endif

#ifndef ROTATINGBOSONSTARSOLUTION_IMPL_HPP_
#define ROTATINGBOSONSTARSOLUTION_IMPL_HPP_

#include <iostream>
#include <vector>
#include <fstream> // Needed for std::ifstream 

RotatingBosonStarSolution::RotatingBosonStarSolution() {}
// inline RotatingBosonStarSolution::RotatingBosonStarSolution(std::string base_path, double BSfreq)
//     : m_base_path(a_base_path), m_BSfreq(a_BSfreq)
// {
// }

void RotatingBosonStarSolution::writeMatrixToFile(const MatDoub& yy, const std::string& filename) {
    std::ofstream outfile(filename);
    
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    int rows = yy.nrows(); // Get the number of rows
    int cols = yy.ncols(); // Get the number of columns

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            outfile << yy[i][j]; // Write the matrix element
            if (j < cols - 1) outfile << " "; // Add a space between elements, avoid trailing space
        }
        outfile << "\n"; // Add a newline after each row
    }

    outfile.close(); // Close the file
    std::cout << "Matrix written to file: " << filename << std::endl;
}

// Read below 
void RotatingBosonStarSolution::main(std::string m_base_path)
{
    std::vector<double> temp_data_x1, temp_data_x2;  // array for theta & array for radius 
    double value_x2, value_x1;
    // std::string base_path = "/Users/macbookpro/Desktop/PHD/Projects/InstabilityRing/";
    std::string filename_radius = "radius.dat";
    std::string filename_theta = "theta.dat";
    std::string filename_amp = "A.dat";
    std::string filename_f = "f.dat";
    std::string filename_l = "l.dat";
    std::string filename_g = "g.dat";
    std::string filename_omega = "omega.dat";
    std::string filename_dthomega = "dthomega.dat"; 
    std::string filename_dromega = "dromega.dat"; 

    std::string full_path_radius = m_base_path + filename_radius;
    std::string full_path_theta = m_base_path + filename_theta;
    std::string full_path_amp = m_base_path + filename_amp;
    std::string full_path_f = m_base_path + filename_f;
    std::string full_path_l = m_base_path + filename_l;
    std::string full_path_g = m_base_path + filename_g;
    std::string full_path_omega = m_base_path + filename_omega;
    std::string full_path_dthomega = m_base_path + filename_dthomega;
    std::string full_path_dromega = m_base_path + filename_dromega;

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
    
    ////////////
    // Theta ///
    ////////////

    std::ifstream file_x2(full_path_theta); 
    if (file_x2.is_open()) {
        // Read each value from the file
        while (file_x2 >> value_x2) {
            temp_data_x2.push_back(value_x2);  // Append to temporary vector
        }
        file_x2.close();
    } else {
        std::cerr << "Unable to open file theta.dat" << std::endl;
    }
    n = temp_data_x2.size();
    x2.resize(n);   
    for (int i = 0; i < n; ++i) {
        x2[i] = temp_data_x2[i];  // Copy from std::vector to VecDoub
    }

    // Resize all matrix variables
    amp.resize(n, m);
    f.resize(n,m);     
    l.resize(n,m);      
    g.resize(n,m);       
    omega.resize(n,m);       
    dthomega.resize(n,m);       
    dromega.resize(n,m);

    ///////////////
    // Amplitude //
    ///////////////

    std::ifstream file_A(full_path_amp);
    if (file_A.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                file_A >> amp[i][j];
            }
        }
        file_A.close();
        } else {
            std::cerr << "Unable to open file A.dat" << std::endl;
        }

    //////////////////
    // f metric var //
    /////////////////

    std::ifstream file_f(full_path_f);
    if (file_f.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                file_f >> f[i][j];
            }
        }
        file_f.close();
        } else {
            std::cerr << "Unable to open file f.dat" << std::endl;
        }

    //////////////////
    // l metric var //
    /////////////////

    std::ifstream file_l(full_path_l);
    if (file_l.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                file_l >> l[i][j];
            }
        }
        file_l.close();
        } else {
            std::cerr << "Unable to open file l.dat" << std::endl;
        }

    //////////////////
    // g metric var //
    /////////////////

    std::ifstream file_g(full_path_g);
    if (file_g.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                file_g >> g[i][j];
            }
        }
        file_g.close();
        } else {
            std::cerr << "Unable to open file g.dat" << std::endl;
        }

    ///////////
    // omega //
    ///////////

    std::ifstream file_omega(full_path_omega);
    if (file_omega.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                file_omega >> omega[i][j];
            }
        }
        file_omega.close();
        } else {
            std::cerr << "Unable to open file omega.dat" << std::endl;
        }

    //////////////
    // dthomega //
    //////////////

    std::ifstream file_dthomega(full_path_dthomega);
    if (file_dthomega.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                file_dthomega >> dthomega[i][j];
            }
        }
        file_dthomega.close();
        } else {
            std::cerr << "Unable to open file dthomega.dat" << std::endl;
        }
        
    //////////////
    // dromega //
    //////////////

    std::ifstream file_dromega(full_path_dromega);
    if (file_dromega.is_open()) {
        // Read each element into the matrix
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                file_dromega >> dromega[i][j];
            }
        }
        file_dromega.close();
        } else {
            std::cerr << "Unable to open file dromega.dat" << std::endl;
        }
}

// 2D interpolation using 1D interpolants
double RotatingBosonStarSolution::interp2D(
    MatDoub Z,              // Z[i][j] values at (x[i], y[j])
    double x_interp,         // Target x
    double y_interp          // Target y
) const
{
    // Step 1: Interpolate along x for each fixed y[j]
    VecDoub intermediate(x2.size());
    for (size_t j = 0; j < x2.size(); ++j) {
        // Extract Z[j][:] (values along x for fixed y[j])
        VecDoub Z_row(Z.ncols());
        for (size_t i = 0; i < x1.size(); ++i) {
            Z_row[i] = Z[j][i];
        }

        Spline_interp cubic_inerpolation(x1,Z_row); 
        intermediate[j] = cubic_inerpolation.interp(x_interp);
        // Interpolate at x_interp
    }

    // Step 2: Interpolate the intermediate values along y
    Spline_interp cubic_inerpolation_new(x2,intermediate); 
    return cubic_inerpolation_new.interp(y_interp);
}

double RotatingBosonStarSolution::get_amp_interp(double r, double theta, int n, int m) const
{
    double interpolated_value; 
    if (x1.size()==0 || x2.size()==0) {
    MayDay::Error("Error: x1 is empty. Make sure to initialize it before calling get_amp_interp().");
    }

    Poly2D_interp bilinear_interpolation(x2,x1,amp,4,4);
    interpolated_value = bilinear_interpolation.interp(theta, r);
    
    // interpolated_value = lagrange_interp_2D(r, theta, x1, m, x2, n, amp);

    return interpolated_value;
}

double RotatingBosonStarSolution::get_f_interp(const double r, const double theta, int n, int m) const
{
    double interpolated_value; 

    Poly2D_interp bilinear_interpolation(x2,x1,f,4,4);
    interpolated_value = bilinear_interpolation.interp(theta, r);

    // interpolated_value = interp2D(f, theta, r);
    
    return interpolated_value;
}

double RotatingBosonStarSolution::get_l_interp(const double r, const double theta, int n, int m) const
{
    double interpolated_value; 

    Poly2D_interp bilinear_interpolation(x2,x1,l,4,4);
    interpolated_value = bilinear_interpolation.interp(theta, r);

    // interpolated_value = interp2D(l, theta, r);
    
    return interpolated_value;
}

double RotatingBosonStarSolution::get_g_interp(const double r, const double theta, int n, int m) const
{
    double interpolated_value; 

    Poly2D_interp bilinear_interpolation(x2,x1,g,4,4);
    interpolated_value = bilinear_interpolation.interp(theta, r);

    // interpolated_value = interp2D(g, theta, r);
    
    return interpolated_value;
}

double RotatingBosonStarSolution::get_omega_interp(const double r, const double theta, int n, int m) const
{
    double interpolated_value; 

    Poly2D_interp bilinear_interpolation(x2,x1,omega,4,4);
    interpolated_value = bilinear_interpolation.interp(theta, r);

    // interpolated_value = interp2D(omega, theta, r);
    
    return interpolated_value;
}

double RotatingBosonStarSolution::get_dthomega_interp(const double r, const double theta, int n, int m) const
{
    double interpolated_value; 
    
    Poly2D_interp bilinear_interpolation(x2,x1,dthomega,4,4);
    interpolated_value = bilinear_interpolation.interp(theta, r);

    // interpolated_value = interp2D(dthomega, theta, r);
    
    return interpolated_value;
}

double RotatingBosonStarSolution::get_dromega_interp(const double r, const double theta, int n, int m) const
{
    double interpolated_value; 

    Poly2D_interp bilinear_interpolation(x2,x1,dromega,4,4);
    interpolated_value = bilinear_interpolation.interp(theta, r);

    // interpolated_value = interp2D(dromega, theta, r);
    
    return interpolated_value;
}

// double RotatingBosonStarSolution::get_BSfrequency() const
// {
//     return BSfreq;
// }

#endif /* ROTATINGBOSONSTARSOLUTION_IMPL_HPP_ */