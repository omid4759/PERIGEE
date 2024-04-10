#ifndef ELEMBC_3D_SLIDING_INTERFACE_HPP
#define ELEMBC_3D_SLIDING_INTERFACE_HPP
// ============================================================================
// ElemBC_3D_SLIDING_INTERFACE.hpp
//
// This is an instantation of ElemBC_3D for the special integral of Nitsche's method
// on the sliding interface.
//
// Since the surface integral involves the first order derivative of basis function,
// we will use the local ien of volume element that has faces on the surface.
// Hence we should write the following data of each boundary (fixed/rotated) into h5 file:
//
// 1. Number of surface elements;
// 2. Face2elem id, i.e. the global id of the volume element where the surface
//    element is attached. We will use it to find local ien of volume element
//    in local assembly.
// 3. Face id of the volume element, we will use the face id and the local ien of
//    volume element to build the quadrature rule and the basis;
//
// To be distinguished from other ElemBC_3D, it uses EBC_Partition_sliding_interface to 
// write h5 file and serves a unique ALocal_interface object.
//
// Author: Xuanming Huang
// Date Created: Apr. 4th 2024
// ============================================================================
#include "ElemBC_3D.hpp"

class ElemBC_3D_sliding_interface : public ElemBC_3D 
{
  public:
    // Input: \para vtkfileList: contains 2 * num_interface_pair vtk files of the interface:
    //              {fixed_0, fixed_1, ..., fixed_n, rotated_0, rotated_1, ..., rotated_n}
    //        \para num_fixed_pt: the number of nodes of the fixed volume
    //              (not only the surface).
    //        \para num_fixed_elem: the number of volume elements of the fixed
    //              volume.
    //              both of the num_fixed_* are used to correct the surface IEN
    //              of the rotated interface.
    //        \para VIEN: IEN of the volume elements.
    //        \para elemtype: element type.
    ElemBC_3D_sliding_interface( const std::vector<std::string> &vtkfileList,
        const int &num_interface_pair_in,
        const int &num_fixed_pt,
        const int &num_fixed_elem,
        const std::vector<double> &vol_ctrlPts,
        const IIEN * const &VIEN,
        const int &elemtype );

    virtual ~ElemBC_3D_sliding_interface();

    virtual int get_num_interface() const
    {return num_interface_pair;}

    virtual int get_fixed_faceID( const int &ii, const int &cell_index ) const
    {return fixed_face_id[ii][cell_index];}

    virtual std::vector<int> get_rotated_faceID(const int &ii) const
    {return rotated_face_id[ii];}

    virtual std::vector<int> get_RL_vien(const int &ii) const
    {return rotated_layer_vien[ii];}

    virtual std::vector<int> get_RLN_GID(const int &ii) const
    {return rotated_layer_global_node[ii];}

    virtual std::vector<double> get_RLN_xyz(const int &ii) const
    {return rotated_layer_pt_xyz[ii];}

  private:
    // the number of [fixed/rotated] interfaces pairs
    const int num_interface_pair;

    // the face id of the volume elements from the fixed interface
    std::vector<std::vector<int>> fixed_face_id;

    // the face id of the volume elements from the rotated interface
    std::vector<std::vector<int>> rotated_face_id;

    // the ien array of the "layer" of rotated volume elements, each of them has a 
    // face on the rotated interface.  
    std::vector<std::vector<int>> rotated_layer_vien;
    
    // the GlobalNodesID of the rotated "layer" volume elements
    std::vector<std::vector<int>> rotated_layer_global_node;

    // the xyz-coordinates of the nodes of the rotated "layer" volume elements
    std::vector<std::vector<double>> rotated_layer_pt_xyz;

    // ------------------------------------------------------------------------
    // Disallow default constructor
    ElemBC_3D_sliding_interface() = delete;
};

#endif
