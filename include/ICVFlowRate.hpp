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
#include <sstream>
#include <string>
#include "Sys_Tools.hpp"

class ICVFlowRate
{
  public:
    ICVFlowRate(){};

    virtual ~ICVFlowRate(){};

    // 0: pulsatile. 1: linear-to-steady. 2: steady.
    virtual int get_inflow_type() const = 0;

    virtual double get_flow_rate( const int &nbc_id, const double &time ) const = 0;

    virtual double get_period( const int &nbc_id ) const
    {
      SYS_T::print_fatal("Error: ICVFlowRate::get_period is not implemented.\n");
      return 0.0;
    }

    virtual int get_num_nbc() const = 0;

    virtual int get_num_bct() const
    {
      SYS_T::print_fatal("Error: ICVFlowRate::get_num_bct is not implemented.\n");
      return 0;
    }

    virtual bool is_bct_id( const int &nbc_id ) const
    {
      SYS_T::print_fatal("Error: ICVFlowRate::is_bct_id is not implemented.\n");
      return false;
    }

    virtual void print_info() const = 0;

  protected:

    // ------------------------------------------------------------------------
    // Generate a filename for inlet face nbc_id as Inlet_xxx_flowrate.txt
    // ------------------------------------------------------------------------
    virtual std::string gen_flowfile_name(const int &nbc_id) const
    {
      std::ostringstream ss;
      ss << "Inlet_";
      if( nbc_id/10 == 0 ) ss << "00";
      else if( nbc_id/100 == 0 ) ss << "0";

      ss << nbc_id << "_precalculated_flowrate.txt";

      return ss.str();
    }
};

#endif
