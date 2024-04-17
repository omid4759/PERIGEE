#include "ALocal_Interface.hpp"

ALocal_Interface::ALocal_Interface( const std::string &fileBaseName, const int &cpu_rank )
: num_itf {0}
{
  const std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/sliding");

  num_itf = h5r -> read_intScalar( gname.c_str(), "num_interface" );
  SYS_T::print_fatal_if(num_itf < 1, "Error, ALocal_Interface: there is no interface in this geometric model.\n");

  num_fixed_ele = h5r -> read_intVector( gname.c_str(), "num_fixed_part_cell" );

  std::string groupbase(gname);
  groupbase.append("/interfaceid_");

  num_rotated_ele.assign(num_itf, 0);
  num_rotated_node.assign(num_itf, 0);

  fixed_vol_ele_id.resize(num_itf);
  fixed_ele_face_id.resize(num_itf);
  rotated_layer_ien.resize(num_itf);
  rotated_layer_face_id.resize(num_itf);
  init_rotated_node_xyz.resize(num_itf);
  rotated_node_id.resize(num_itf);

  for(int ii=0; ii<num_itf; ++ii)
  {
    if( num_fixed_ele[ii] > 0 )
    {
      std::string subgroup_name(groupbase);
      subgroup_name.append( std::to_string(ii) );

      fixed_vol_ele_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "part_fixed_cell_id" );

      fixed_ele_face_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "fixed_cell_face_id" );

      rotated_layer_ien[ii] = h5r -> read_intVector( subgroup_name.c_str(), "rotated_cell_ien" );

      rotated_layer_face_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "rotated_cell_face_id" );

      init_rotated_node_xyz[ii] = h5r -> read_doubleVector( subgroup_name.c_str(), "rotated_node_xyz" );

      rotated_node_id[ii] = h5r -> read_intVector( subgroup_name.c_str(), "rotated_node_map" );
    }
    else
    {
      fixed_vol_ele_id[ii].clear();
      fixed_ele_face_id[ii].clear();
      rotated_layer_ien[ii].clear();
      rotated_layer_face_id[ii].clear();
      init_rotated_node_xyz[ii].clear();
      rotated_node_id[ii].clear();
    }
  }

  const std::string mesh_info("Global_Mesh_Info");
  nLocBas = h5r -> read_intScalar(mesh_info.c_str(), "nLocBas");

  delete h5r; H5Fclose( file_id );
}

void ALocal_Interface::print_info() const
{
  SYS_T::commPrint("Interfaces: %d\n", num_itf);
}

Vector_3 ALocal_Interface::get_curr_xyz(const int &ii, const int &node, const double &tt) const
{
  Vector_3 xyz (get_init_rotated_node_xyz(ii, 3 * node),
                get_init_rotated_node_xyz(ii, 3 * node + 1),
                get_init_rotated_node_xyz(ii, 3 * node + 2));

  // rotation around z-axis
  const double angular_velo = MATH_T::PI / 60;  // (rad/s)

  const double rr = std::sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2));

  double angle = MATH_T::get_angle_2d(xyz(1), xyz(2));

  angle += angular_velo * tt;

  xyz(1) = std::cos(angle) * rr;
  xyz(2) = std::sin(angle) * rr;

  return xyz;
}

void ALocal_Interface::get_ele_ctrlPts(const int &ii, const int &ee, const double &tt,
  double * const volctrl_x,  double * const volctrl_y,  double * const volctrl_z) const
{
  for(int nn{0}; nn < nLocBas; ++nn)
  {
    int node = get_rotated_layer_ien(ii, nLocBas * ee + nn);
    Vector_3 cuur_node_xyz = get_curr_xyz(ii, node, tt);

    volctrl_x[nn] = cuur_node_xyz(0);
    volctrl_y[nn] = cuur_node_xyz(1);
    volctrl_z[nn] = cuur_node_xyz(2);
  }
}

// EOF