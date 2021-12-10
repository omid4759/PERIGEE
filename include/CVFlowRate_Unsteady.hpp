#ifndef CVFLOWRATE_UNSTEADY_HPP
#define CVFLOWRATE_UNSTEADY_HPP
// ==================================================================
// CVFlowRate_Unsteady.hpp
//
// Unsteady flow for inflow condition.
// 
// The period value should be compatible with the original data file.
//
// I use time / period to determine the past full cycles. Then
// time - num_of_past_period * period gives the time in the local
// cycle, which can be used to evaluate the fourier wave.
//
// Output flow rate is defined as follows
//   a_0 + sum_i a_i cos(i w t) + b_i sin(i w t) 
//       for i = 1 : num_of_mode
//           t = t - [t/period] x period
//   
// Author: Ju Liu
// Date Created: Sept. 23 2017
// ==================================================================
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "ICVFlowRate.hpp"

class CVFlowRate_Unsteady : public ICVFlowRate
{
  public:
    CVFlowRate_Unsteady( const std::string &filename );

    virtual ~CVFlowRate_Unsteady();

    virtual double get_flow_rate( const int &nbc_id, const double &time ) const;

    virtual double get_period( const int &nbc_id ) const { return period[nbc_id]; };

    virtual int get_num_nbc() const { return num_nbc; }

    virtual int get_num_bct() const { return bct_id.size(); }

    virtual bool is_bct_id( const int &nbc_id ) const
    { return VEC_T::is_invec(bct_id, nbc_id); }

    virtual void print_info() const;

  private:
    int num_nbc;

    // nbc_id's for which to assign velocity profiles
    // length num_bct
    std::vector<int> bct_id;

    std::vector< std::vector<double> > coef_a, coef_b;

    std::vector<int> num_of_mode;

    std::vector<double> w, period;
};

#endif
