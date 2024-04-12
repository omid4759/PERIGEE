#include "EBC_Partition_sliding_interface.hpp"

EBC_Partition_sliding_interface::EBC_Partition_sliding_interface(const IPart * const &part,
    const Map_Node_Index * const &mnindex, const ElemBC * const &ebc)
: EBC_Partition(part, mnindex, ebc), num_pair{ebc->get_num_interface()}, num_fixed_part_ele{std::vector<int> {}},
  fixed_part_vol_ele_id{std::vector<std::vector<int>> {}}, fixed_ele_face_id{std::vector<std::vector<int>> {}},
  rotated_layer_ien{std::vector<std::vector<int>> {}}, rotated_ele_face_id{std::vector<std::vector<int>> {}},
  rotated_layer_global_node{std::vector<std::vector<int>> {}}, rotated_layer_pt_xyz{std::vector<std::vector<double>> {}}
{
  num_fixed_part_ele.resize(num_pair);
  fixed_part_vol_ele_id.resize(num_pair);
  fixed_ele_face_id.resize(num_pair);
  rotated_layer_ien.resize(num_pair);
  rotated_ele_face_id.resize(num_pair);
  rotated_layer_global_node.resize(num_pair);
  rotated_layer_pt_xyz.resize(num_pair);

  for(int ii=0; ii<num_pair; ++ii)
  {
    num_fixed_part_ele[ii] = num_local_cell[ii];

    fixed_part_vol_ele_id[ii] = std::vector<int> {};
    fixed_ele_face_id[ii] = std::vector<int> {};
    PERIGEE_OMP_PARALLEL
    {
      std::vector<int> temp_part_vol_ele_id {};
      std::vector<int> temp_ele_face_id {};

      PERIGEE_OMP_FOR
      for(int ee=0; ee < ebc->get_num_cell(ii); ++ee)
      {
        const int global_vol_ele_id = ebc -> get_global_cell(ii, ee);
        const int loc_id = part -> get_elemLocIndex(global_vol_ele_id);
        if( loc_id != -1 )
        {
          temp_part_vol_ele_id.push_back( loc_id );
          temp_ele_face_id.push_back( ebc -> get_fixed_faceID(ii, ee) );
        }
      }

      PERIGEE_OMP_CRITICAL
      {
        VEC_T::insert_end(fixed_part_vol_ele_id[ii], temp_part_vol_ele_id);
        VEC_T::insert_end(fixed_ele_face_id[ii], temp_ele_face_id);
      }
    }

    rotated_layer_ien[ii] = std::vector<int> {};
    rotated_ele_face_id[ii] = std::vector<int> {};
    rotated_layer_global_node[ii] = std::vector<int> {};
    rotated_layer_pt_xyz[ii] = std::vector<double> {};

    // write the rotated layer's info only when this part has fixed interface element
    if(num_local_cell[ii] > 0)
    {
      rotated_layer_ien[ii] = ebc -> get_RL_vien(ii);

      rotated_ele_face_id[ii] = ebc -> get_rotated_faceID(ii);
      
      rotated_layer_global_node[ii] = ebc -> get_RLN_GID(ii);

      // convert the GlobalNodeID to new(mapped) global node id
      PERIGEE_OMP_PARALLEL_FOR
      for(int jj=0; jj<VEC_T::get_size(rotated_layer_global_node[ii]); ++jj)
      {
        const int new_gid = mnindex->get_old2new(rotated_layer_global_node[ii][jj]);
        rotated_layer_global_node[ii][jj] = new_gid;
      }

      rotated_layer_pt_xyz[ii] = ebc -> get_RLN_xyz(ii);
    }
    else
      ; // do nothing
  }
}

EBC_Partition_sliding_interface::~EBC_Partition_sliding_interface()
{
  for(int ii=0; ii<num_pair; ++ii)
  {
    fixed_part_vol_ele_id[ii].clear();
    fixed_ele_face_id[ii].clear();
    rotated_layer_ien[ii].clear();
    rotated_ele_face_id[ii].clear();
    rotated_layer_global_node[ii].clear();
    rotated_layer_pt_xyz[ii].clear();
  }
  fixed_part_vol_ele_id.clear();
  fixed_ele_face_id.clear();
  rotated_layer_ien.clear();
  rotated_ele_face_id.clear();
  rotated_layer_global_node.clear();
  rotated_layer_pt_xyz.clear();
}

void EBC_Partition_sliding_interface::write_hdf5(const std::string &FileName) const
{
  const std::string fName = SYS_T::gen_partfile_name( FileName, cpu_rank );

  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t g_id = H5Gcreate(file_id, "/sliding", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  h5w -> write_intScalar( g_id, "num_interface", num_pair );

  h5w -> write_intVector( g_id, "num_part_fixed_cell", num_fixed_part_ele );

  const std::string groupbase("interfaceid_");

  for(int ii=0; ii<num_pair; ++ii)
  {
    std::string subgroup_name(groupbase);
    subgroup_name.append( std::to_string(ii) );

    hid_t group_id = H5Gcreate(g_id, subgroup_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if(num_local_cell[ii] > 0)
    {
      h5w -> write_intVector( group_id, "part_fixed_cell_id", fixed_part_vol_ele_id[ii] );

      h5w -> write_intVector( group_id, "fixed_cell_face_id", fixed_ele_face_id[ii] );

      h5w -> write_intScalar( group_id, "num_rotated_cell", VEC_T::get_size(rotated_ele_face_id[ii]) );

      h5w -> write_intScalar( group_id, "num_rotated_node", VEC_T::get_size(rotated_layer_global_node[ii]) );

      h5w -> write_intVector( group_id, "rotated_cell_ien", rotated_layer_ien[ii] );

      h5w -> write_intVector( group_id, "rotated_cell_face_id", rotated_ele_face_id[ii] );

      h5w -> write_intVector( group_id, "rotated_node_map", rotated_layer_global_node[ii] );

      h5w -> write_doubleVector( group_id, "rotated_node_xyz", rotated_layer_pt_xyz[ii] );
    }
    else
      ; // stop writing if num_part_fixed_cell = 0

    H5Gclose( group_id );
  }
  delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
}

// EOF
