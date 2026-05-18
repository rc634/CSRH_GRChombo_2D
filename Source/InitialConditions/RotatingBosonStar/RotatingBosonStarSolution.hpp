#ifndef ROTATINGBOSONSTARSOLUTION_HPP_
#define ROTATINGBOSONSTARSOLUTION_HPP_

#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "nr3.h"
#include "interp_1d.hpp"
#include "interp_linear.hpp"
#include "interp_2d.hpp"

class RotatingBosonStarSolution
{
  // protected:        
  //   double m_BSfreq;
  //   std::string m_base_path;
  public:
    RotatingBosonStarSolution();
    // RotatingBosonStarSolution(std::string a_base_path, double a_BSfreq);
    int m, n;             // size of matrix used 
    VecDoub x1, x2;       // array for theta & array for radius 
    MatDoub amp;            // scalar field modulus
    MatDoub f;           // f metric function
    MatDoub l;          // l metric function
    MatDoub g;         // g metric function
    MatDoub omega;        // omega metric function
    MatDoub dthomega;       // theta derivative of omega
    MatDoub dromega;       // r derivative of omega

    double get_amp_interp(const double r, const double theta, int n, int m) const;
    double get_f_interp(const double r, const double theta, int n, int m) const;
    double get_l_interp(const double r, const double theta, int n, int m) const;
    double get_g_interp(const double r, const double theta, int n, int m) const;
    double get_omega_interp(const double r, const double theta, int n, int m) const;
    double get_dthomega_interp(const double r, const double theta, int n, int m) const;
    double get_dromega_interp(const double r, const double theta, int n, int m) const;
    double interp2D(MatDoub Z, double x_interp, double y_interp) const;
    // double get_BSfrequency() const;
    void writeMatrixToFile(const MatDoub& yy, const std::string& filename);
    void main(std::string m_base_path);
}; 

#include "RotatingBosonStarSolution.impl.hpp"

#endif /* ROTATINGBOSONSTARSOLUTION_HPP_ */