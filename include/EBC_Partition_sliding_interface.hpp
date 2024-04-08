#ifndef EBC_PARTITION_SLIDING_INTERFACE_HPP
#define EBC_PARTITION_SLIDING_INTERFACE_HPP
// ==================================================================
// EBC_Partition_sliding_interface.hpp
//
// Element boundary condition partition for Nitsche's method applied
// on the interface problem.
// 
// Author: Xuanming Huang
// Date: Apr. 4th 2023
// ==================================================================
#include "EBC_Partition.hpp"

class EBC_Partition_sliding_interface : public EBC_Partition
{
  public:
    // The input ElemBC should be ElemBC_3D_wall_turbulence
    EBC_Partition_sliding_interface( const IPart * const &part,
        const Map_Node_Index * const &mnidex,
        const ElemBC * const &ebc );

    virtual ~EBC_Partition_sliding_interface();

    // write the data to hdf5 file in folder /interface
    virtual void write_hdf5( const std::string &FileName ) const;

  protected:
    // stores the local volume element id of the fixed interface in this part
    std::vector<int> fixed_part_vol_ele_id;

    // stores the face id of the volume element of the fixed interface in this part
    std::vector<int> fixed_ele_face_id;

    // stores ien of all the volume element of the rotated layer
    std::vector<int> rotated_layer_ien;

    // stores the face id of all the volume element of the rotated layer
    std::vector<int> rotated_ele_face_id;

    // stores the rotated layer nodes indices, converted by get_old2new()
    std::vector<int> rotated_layer_nodes;

    // stores the rotated layer nodes' coordinates
    std::vector<double> rotated_layer_nodes_xyz;
};

#endif
