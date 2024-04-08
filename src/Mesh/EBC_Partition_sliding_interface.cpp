#include "EBC_Partition_sliding_interface.hpp"

EBC_Partition_sliding_interface::EBC_Partition_sliding_interface(const IPart * const &part,
    const Map_Node_Index * const &mnindex, const ElemBC * const &ebc)
: EBC_Partition(part, mnindex, ebc),
  fixed_part_vol_ele_id{std::vector<int> {}}, fixed_ele_face_id{std::vector<int> {}},
  rotated_layer_ien{std::vector<int> {}}, rotated_ele_face_id{std::vector<int> {}},
  rotated_layer_nodes{std::vector<int> {}}, rotated_layer_nodes_xyz{std::vector<double> {}}
{
  PERIGEE_OMP_PARALLEL
  {
    std::vector<int> temp_part_vol_ele_id {};
    std::vector<int> temp_ele_face_id {};

    PERIGEE_OMP_FOR
    for(int ee=0; ee < ebc->get_num_cell(0); ++ee)
    {
      const int global_vol_ele_id = ebc -> get_global_cell(0, ee);
      const int loc_id = part -> get_elemLocIndex(global_vol_ele_id);
      if( loc_id != -1 )
      {
        temp_part_vol_ele_id.push_back( loc_id );
        temp_ele_face_id.push_back( ebc -> get_ifaceID(0, ee) );
      }
    }

    PERIGEE_OMP_CRITICAL
    {
      VEC_T::insert_end(fixed_part_vol_ele_id, temp_part_vol_ele_id);
      VEC_T::insert_end(fixed_ele_face_id, temp_ele_face_id);
    }
  }
  
  // write the rotated layer's info only when this part has fixed interface element
  if(ebc->get_num_cell(0) > 0)
  {
    rotated_layer_ien = ebc -> get_vien_RL();

    rotated_ele_face_id.resize(ebc->get_num_cell(1));
    PERIGEE_OMP_PARALLEL_FOR
    for(int ee=0; ee < ebc->get_num_cell(1); ++ee)
      rotated_ele_face_id[ee] = ebc -> get_ifaceID(1, ee);
    
    rotated_layer_nodes = ebc -> get_RLN_GID();
    PERIGEE_OMP_PARALLEL_FOR
    for(int ii=0; ii < VEC_T::get_size(rotated_layer_nodes); ++ii)
    {
      const int new_gid = mnindex->get_old2new(rotated_layer_nodes[ii]);
      rotated_layer_nodes[ii] = new_gid;
    }

    rotated_layer_nodes_xyz = ebc -> get_RLN_xyz();
  }
  else
    ; // do nothing
}

EBC_Partition_sliding_interface::~EBC_Partition_sliding_interface()
{
  fixed_part_vol_ele_id.clear();
  fixed_ele_face_id.clear();
  rotated_layer_ien.clear();
  rotated_ele_face_id.clear();
  rotated_layer_nodes.clear();
  rotated_layer_nodes_xyz.clear();
}

void EBC_Partition_sliding_interface::write_hdf5(const std::string &FileName) const
{
  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, "/sliding", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  h5w -> write_intScalar( g_id, "num_part_fixed_cell", num_local_cell[0] );

  h5w -> write_intVector( g_id, "part_fixed_cell_id", fixed_part_vol_ele_id );

  h5w -> write_intVector( g_id, "fixed_cell_face_id", fixed_ele_face_id );

  if(num_local_cell[0] > 0)
  {
    h5w -> write_intScalar( g_id, "num_all_rotated_cell", VEC_T::get_size(rotated_ele_face_id) );

    h5w -> write_intScalar( g_id, "num_all_rotated_node", VEC_T::get_size(rotated_layer_nodes) );

    h5w -> write_intVector( g_id, "all_rotated_cell_ien", rotated_layer_ien );

    h5w -> write_intVector( g_id, "rotated_cell_face_id", rotated_ele_face_id );

    h5w -> write_intVector( g_id, "rotated_node_map", rotated_layer_nodes );

    h5w -> write_doubleVector( g_id, "rotated_node_xyz", rotated_layer_nodes_xyz );
  }
  else
    ;   // stop writing if num_part_fixed_cell = 0

  delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
}

// EOF
