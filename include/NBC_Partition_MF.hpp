#ifndef NBC_PARTITION_MF_HPP
#define NBC_PARTITION_MF_HPP
// ============================================================================
// NBC_Partition_MF.hpp
//
// Nodal boundary condition implementation for Multi-Field problems.
// 
// For Multi-Field problems, we need to document the assigned row (or column)
// index for the (node, dof) pair.
//
// Date: Dec. 17 2021
// Author: Ju Liu
// ============================================================================
#include "NBC_Partition.hpp"

class NBC_Partition_MF : public NBC_Partition
{
  public:
    NBC_Partition_MF( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const std::vector<INodalBC *> &nbc_list,
        const std::vector< std::vector<int> > &grid2id );

    virtual ~NBC_Partition_MF();

  protected:
    // LID mapped to the MF value.
    std::vector<int> LID_MF;

    // LDN, LPSN, LPMN, LocalMaster, LocalMasterSlave mapped to MF value
    // The mapping is from the grid2id mapper, that maps the (node, dof) value
    // to the actual matrix problem's row/col location.
    std::vector<int> LDN_MF, LPSN_MF, LPMN_MF, LocalMaster_MF, LocalMasterSlave_MF;
};

#endif
