#ifndef NBC_PARTITION_3D_INFLOW_HPP
#define NBC_PARTITION_3D_INFLOW_HPP
// ==================================================================
// NBC_Partition_3D_inflow.hpp
//
// Inflow Nodal Boundary condition partition implementation for 3D
// meshes. 
//
// This NBC partition code is specifically designed for the class of
// NodalBC_3D_inflow, which contains two additional information:
// inflow surface area and the inflow surface outward normal vector.
// 
// The data recorded in the HDF5 file by this class will be loaded
// in the ALocal_Inflow_NodalBC class in the analysis code.
//
// Date crated: Aug. 9 2017
// Author: Ju Liu
// ==================================================================
#include "NBC_Partition_3D.hpp"

class NBC_Partition_3D_inflow : public NBC_Partition_3D
{
  public:
    NBC_Partition_3D_inflow(const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const INodalBC * const &nbc );

    virtual ~NBC_Partition_3D_inflow();

    virtual void write_hdf5( const std::string &FileName ) const;

  private:
    // active area & full area. Length num_nbc
    std::vector<double> actarea, facearea;

    // unit outward normal vector on this surface. Length num_nbc
    std::vector<Vector_3> outvec;

    // number of boundary points of this surface. Length num_nbc
    std::vector<int> num_out_bc_pts;

    // spatial coordinates of the centroid. Length num_nbc
    std::vector<Vector_3> centroid;

    // coordinates of the boundary points
    // num_nbc x num_out_bc_pts[ii]
    std::vector< std::vector<double> > outline_pts;

    // number of local node / cell, element number of nodes. Length num_nbc
    std::vector<int> num_local_node, num_local_cell, cell_nLocBas;

    // local nodes' coordinates
    // num_nbc times (3 x num_local_node[ii])
    std::vector< std::vector<double> > local_pt_xyz;

    // local cell IEN array
    // num_nbc times (cell_nLocBas[ii] x num_local_cell[ii])
    std::vector< std::vector<int> > local_tri_ien;

    // local node's global index
    // num_nbc x num_local_node[ii]
    std::vector< std::vector<int> > local_global_node;

    // local node's position in the local_to_global_array
    // num_nbc x num_local_node[ii]
    std::vector< std::vector<int> > local_node_pos;

    // local cell's global index
    // num_nbc x num_local_cell[ii]
    std::vector< std::vector<int> > local_global_cell;
  
    // ------------------------------------------------------------------------
    // This function is NOT allowed for inflow nbc 
    virtual void write_hdf5( const std::string &FileName,
        const std::string &GroupName ) const;
};

#endif
