#ifndef SINGLEBOSONSTAR_HPP_
#define SINGLEBOSONSTAR_HPP_

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

class SingleBosonStar
{
  public:
    SingleBosonStar();

    int m;    
    VecDoub x1;        
    VecDoub amp;         
    VecDoub f;          
    VecDoub Phi;          

    double get_amp_interp(const double r, int m) const;
    double get_f_interp(const double r, int m) const;
    double get_lapse_interp(const double r, int m) const;
    void main(std::string m_base_path);
}; 

#include "SingleBosonStar.impl.hpp"

#endif /* SINGLEBOSONSTAR_HPP_ */
