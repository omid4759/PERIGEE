#ifndef ICVFLOWRATE_HPP
#define ICVFLOWRATE_HPP
// ==================================================================
// ICVFlowRate.hpp
//
// Interface file for the cardiovascular inflow flow rate function.
//
// This function shall have two instantiations: steady case and 
// unsteady case.
//
// Author: Ju Liu
// Date created: Aug. 6 2017
// ==================================================================

#include "Sys_Tools.hpp"

class ICVFlowRate
{
  public:
    ICVFlowRate(){};

    virtual ~ICVFlowRate(){};

    virtual double get_flow_rate( const double &time ) const = 0;

    virtual int get_velo_profile_type() const
    { SYS_T::commPrint("Warning: get_velo_profile_type is not implemented. \n"); return -1; }

    virtual void get_fourier_coeff( std::vector<double>& a_n, std::vector<double>& b_n ) const
    { SYS_T::commPrint("Warning: get_fourier_coeff is not implemented. \n"); }

    virtual void print_info() const = 0;
};

#endif
