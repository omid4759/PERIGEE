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
    // Input: \para vtkfileList: contains two vtk files of the interface.
    //              The first surface is attached to the fixed volume, 
    //              and the second one is attached the the rotated volume.
    //        \para num_fixed_pt: the number of nodes of the fixed volume
    //              (not only the surface).
    //        \para num_fixed_elem: the number of volume elements of the fixed
    //              volume.
    //              both of the num_fixed_* are used to correct the surface IEN
    //              of the rotated interface.
    //        \para VIEN: IEN of the volume elements.
    //        \para elemtype: element type.
    ElemBC_3D_sliding_interface( const std::vector<std::string> &vtkfileList,
        const int &num_fixed_pt,
        const int &num_fixed_elem,
        const std::vector<double> &vol_ctrlPts,
        const IIEN * const &VIEN,
        const int &elemtype );

    virtual ~ElemBC_3D_sliding_interface();

    virtual int get_ifaceID( const int &f_or_r, const int &cell_index ) const
    {return face_id[f_or_r][cell_index];}

    virtual std::vector<int> get_vien_RL() const
    {return vien_rotated_layer;}

    virtual std::vector<int> get_RLN_GID() const
    {return layer_nodes_GID;}

    virtual std::vector<double> get_RLN_xyz() const
    {return layer_nodes_xyz;}

  private:
    // the face id of the volume element
    // face_id[0][*]: the fixed interface
    // face_id[1][*]: the rotated interface
    std::vector<std::vector<int>> face_id;

    // the ien array of the "layer" of rotated volume elements, each of them has a 
    // face on the rotated interface.  
    std::vector<int> vien_rotated_layer;
    
    // the GlobalNodesID of the rotated "layer" volume elements
    std::vector<int> layer_nodes_GID;

    // the xyz-coordinates of the nodes of the rotated "layer" volume elements
    std::vector<double> layer_nodes_xyz;

    // ------------------------------------------------------------------------
    // Disallow default constructor
    ElemBC_3D_sliding_interface() = delete;
};

#endif
