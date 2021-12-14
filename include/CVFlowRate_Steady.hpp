#ifndef CVFLOWRATE_STEADY_HPP
#define CVFLOWRATE_STEADY_HPP
// ============================================================================
// CVFlowRate_Steady.hpp
//
// Steady flow for inflow condition
//
// Author: Ju Liu
// Date Created: May 1 2021
// ============================================================================
#include "Vec_Tools.hpp"
#include "Math_Tools.hpp"
#include "ICVFlowRate.hpp"

class CVFlowRate_Steady : public ICVFlowRate
{
  public:
    CVFlowRate_Steady( const std::string &filename );

    CVFlowRate_Steady( const int &input_num_nbc, const double &in_flowrate );

    virtual ~CVFlowRate_Steady();

    virtual int get_inflow_type() const { return 1; }

    virtual double get_flow_rate(const int &nbc_id, const double &time) const;

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

    std::vector<double> flowrate;
};

#endif
