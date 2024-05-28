#include "ElemBC_3D_sliding_interface.hpp"

ElemBC_3D_sliding_interface::ElemBC_3D_sliding_interface( 
    const std::vector<std::string> &vtkfileList,
    const int &num_interface_pair_in,
    const int &num_fixed_pt, const int &num_fixed_elem,
    const std::vector<double> &vol_ctrlPts, const IIEN * const &VIEN, const int &elemtype )
: ElemBC_3D ( vtkfileList, elemtype ), num_interface_pair{num_interface_pair_in},
  fixed_face_id{std::vector<std::vector<int>> {}}, fixed_part_tag{std::vector<std::vector<int>> {}},
  fixed_layer_vien{std::vector<std::vector<int>> {}}, fixed_layer_global_node{std::vector<std::vector<int>> {}},
  fixed_layer_pt_xyz{std::vector<std::vector<double>> {}}, fixed_interval_tag{std::vector<std::vector<int>> {}},
  rotated_face_id{std::vector<std::vector<int>> {}},
  rotated_layer_vien{std::vector<std::vector<int>> {}}, rotated_layer_global_node{std::vector<std::vector<int>> {}},
  rotated_layer_pt_xyz{std::vector<std::vector<double>> {}}, rotated_interval_tag{std::vector<std::vector<int>> {}}
{
  SYS_T::print_fatal_if(VEC_T::get_size(vtkfileList) != 2 * num_interface_pair,
    "Error, ElemBC_3D_sliding_interface: The number of interface file is wrong!\n");
  
  fixed_face_id.resize(num_interface_pair);
  fixed_part_tag.resize(num_interface_pair);
  fixed_layer_vien.resize(num_interface_pair);
  rotated_face_id.resize(num_interface_pair);
  rotated_layer_vien.resize(num_interface_pair);
  for(int ii=0; ii<num_interface_pair; ++ii)
  {
    fixed_face_id[ii].resize(num_cell[ii]);
    fixed_layer_vien[ii].resize(VIEN->get_nLocBas() * num_cell[ii]);

    rotated_face_id[ii].resize(num_cell[ii + num_interface_pair]);
    rotated_layer_vien[ii].resize(VIEN->get_nLocBas() * num_cell[ii + num_interface_pair]);
  }

  fixed_layer_global_node.resize(num_interface_pair);
  fixed_layer_pt_xyz.resize(num_interface_pair);
  fixed_interval_tag.resize(num_interface_pair);
  rotated_layer_global_node.resize(num_interface_pair);
  rotated_layer_pt_xyz.resize(num_interface_pair);
  rotated_interval_tag.resize(num_interface_pair);

  const std::string filename_base = "epart_";
  const std::string filename_tail = "_itf.h5";

  if(elem_type == 501 || elem_type == 502)
  {
    TET_T::Tet4 * tetcell = new TET_T::Tet4();

    for(int ii=0; ii<num_interface_pair; ++ii)
    {
      const std::string filename = SYS_T::gen_capfile_name(filename_base, ii, filename_tail);
      fixed_part_tag[ii] = HDF5_T::read_intVector( filename.c_str(), "/", "part" );

      // Find the face id of the fixed interfaces
      for(int ee=0; ee<num_cell[ii]; ++ee)
      {
        const int node_t[3] { get_ien(ii, ee, 0), get_ien(ii, ee, 1), get_ien(ii, ee, 2) };

        const int node_t_gi[3] { get_global_node(ii, node_t[0]),
                                 get_global_node(ii, node_t[1]),
                                 get_global_node(ii, node_t[2]) };

        const int cell_gi = get_global_cell(ii, ee);

        const int tet_n[4] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3) };
        
        tetcell->reset(tet_n[0], tet_n[1], tet_n[2], tet_n[3]);

        fixed_face_id[ii][ee] = tetcell->get_face_id(node_t_gi[0], node_t_gi[1], node_t_gi[2]);

        // Construct the "layer" vien of the fixed interfaces
        for(int kk=0; kk < VIEN->get_nLocBas(); ++kk)
          fixed_layer_vien[ii][ee * VIEN->get_nLocBas() + kk] = VIEN->get_IEN(cell_gi, kk);
      }

      // Find the face id of the rotated interfaces
      const int jj = ii + num_interface_pair;
      for(int ee=0; ee<num_cell[jj]; ++ee)
      {
        const int node_t[3] { get_ien(jj, ee, 0), get_ien(jj, ee, 1), get_ien(jj, ee, 2) };

        const int node_t_gi[3] { get_global_node(jj, node_t[0]) + num_fixed_pt,
                                 get_global_node(jj, node_t[1]) + num_fixed_pt,
                                 get_global_node(jj, node_t[2]) + num_fixed_pt};

        global_cell[jj][ee] += num_fixed_elem;
        const int cell_gi = get_global_cell(jj, ee);

        const int tet_n[4] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3) };
        
        tetcell->reset(tet_n[0], tet_n[1], tet_n[2], tet_n[3]);

        rotated_face_id[ii][ee] = tetcell->get_face_id(node_t_gi[0], node_t_gi[1], node_t_gi[2]);

        // Construct the "layer" vien of the rotated interfaces
        for(int kk=0; kk < VIEN->get_nLocBas(); ++kk)
          rotated_layer_vien[ii][ee * VIEN->get_nLocBas() + kk] = VIEN->get_IEN(cell_gi, kk);
      }
    }
    delete tetcell;

    std::vector<double> intervals_0 {-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
    std::vector<double> intervals_12 {0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.30};

    Interface_Element_Group IEG_fixed_0(vtkfileList[0], get_cell_nLocBas(0), 0, intervals_0);
    fixed_interval_tag[0] = IEG_fixed_0.get_tag();
    Interface_Element_Group IEG_rotated_0(vtkfileList[3], get_cell_nLocBas(3), 0, intervals_0);
    rotated_interval_tag[0] = IEG_rotated_0.get_tag();

    Interface_Element_Group IEG_fixed_1(vtkfileList[1], get_cell_nLocBas(1), Vector_3(0.5, 0.0, 0.0), intervals_12);
    fixed_interval_tag[1] = IEG_fixed_1.get_tag();

    Interface_Element_Group IEG_rotated_1(vtkfileList[4], get_cell_nLocBas(4), Vector_3(0.5, 0.0, 0.0), intervals_12);
    rotated_interval_tag[1] = IEG_rotated_1.get_tag();


    Interface_Element_Group IEG_fixed_2(vtkfileList[2], get_cell_nLocBas(2), Vector_3(-0.5, 0.0, 0.0), intervals_12);
    fixed_interval_tag[2] = IEG_fixed_2.get_tag();

    Interface_Element_Group IEG_rotated_2(vtkfileList[5], get_cell_nLocBas(5), Vector_3(-0.5, 0.0, 0.0), intervals_12);
    rotated_interval_tag[2] = IEG_rotated_2.get_tag();
  }
  else if(elem_type == 601 || elem_type == 602)
  {
    HEX_T::Hex8 * hexcell = new HEX_T::Hex8();

    for(int ii=0; ii<num_interface_pair; ++ii)
    {
      const std::string filename = SYS_T::gen_capfile_name(filename_base, ii, filename_tail);
      fixed_part_tag[ii] = HDF5_T::read_intVector( filename.c_str(), "/", "part" );

      // Find the face id of the fixed interfaces
      for(int ee=0; ee<num_cell[ii]; ++ee)
      {
        const int node_q[4] { get_ien(ii, ee, 0), get_ien(ii, ee, 1),
                              get_ien(ii, ee, 2), get_ien(ii, ee, 3) };
        
        const int node_q_gi[4] { get_global_node(ii, node_q[0]),
                                 get_global_node(ii, node_q[1]),
                                 get_global_node(ii, node_q[2]),
                                 get_global_node(ii, node_q[3]) };
        
        const int cell_gi = get_global_cell(ii, ee);

        const int hex_n[8] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3),
                             VIEN->get_IEN(cell_gi, 4), VIEN->get_IEN(cell_gi, 5),
                             VIEN->get_IEN(cell_gi, 6), VIEN->get_IEN(cell_gi, 7) };
        
        hexcell->reset(hex_n[0], hex_n[1], hex_n[2], hex_n[3],
                       hex_n[4], hex_n[5], hex_n[6], hex_n[7]);

        fixed_face_id[ii][ee] = hexcell->get_face_id(node_q_gi[0], node_q_gi[1], node_q_gi[2], node_q_gi[3]);

        // Construct the "layer" vien of the fixed interfaces
        for(int kk=0; kk < VIEN->get_nLocBas(); ++kk)
          fixed_layer_vien[ii][ee * VIEN->get_nLocBas() + kk] = VIEN->get_IEN(cell_gi, kk);
      }

      // Find the face id of the rotated interfaces
      const int jj = ii + num_interface_pair;
      for(int ee=0; ee<num_cell[jj]; ++ee)
      {
        const int node_q[4] { get_ien(jj, ee, 0), get_ien(jj, ee, 1),
                              get_ien(jj, ee, 2), get_ien(jj, ee, 3) };
        
        const int node_q_gi[4] { get_global_node(jj, node_q[0]) + num_fixed_pt,
                                 get_global_node(jj, node_q[1]) + num_fixed_pt,
                                 get_global_node(jj, node_q[2]) + num_fixed_pt,
                                 get_global_node(jj, node_q[3]) + num_fixed_pt};
        
        global_cell[jj][ee] += num_fixed_elem;
        const int cell_gi = get_global_cell(jj, ee);

        const int hex_n[8] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3),
                             VIEN->get_IEN(cell_gi, 4), VIEN->get_IEN(cell_gi, 5),
                             VIEN->get_IEN(cell_gi, 6), VIEN->get_IEN(cell_gi, 7) };
        
        hexcell->reset(hex_n[0], hex_n[1], hex_n[2], hex_n[3],
                       hex_n[4], hex_n[5], hex_n[6], hex_n[7]);

        rotated_face_id[ii][ee] = hexcell->get_face_id(node_q_gi[0], node_q_gi[1], node_q_gi[2], node_q_gi[3]);

        // Construct the "layer" vien of the rotated interfaces
        for(int kk=0; kk<VIEN->get_nLocBas(); ++kk)
          rotated_layer_vien[ii][ee * VIEN->get_nLocBas() + kk] = VIEN->get_IEN(cell_gi, kk);
      }
    }
    delete hexcell;
  }
  else
    SYS_T::print_fatal("Error: ElemBC_3D_sliding_interface, unknown element type.\n");
  
  for(int ii=0; ii<num_interface_pair; ++ii)
  {
    fixed_layer_global_node[ii] = fixed_layer_vien[ii];
    VEC_T::sort_unique_resize(fixed_layer_global_node[ii]);
    const int num_fixed_layer_node = VEC_T::get_size(fixed_layer_global_node[ii]);

    // Convert the GlobalNodeID in fixed_layer_vien to local indices in this ebc
    PERIGEE_OMP_PARALLEL_FOR
    for(int &nodeid : fixed_layer_vien[ii])
    {
      const int local_id = VEC_T::get_pos(fixed_layer_global_node[ii], nodeid);
      nodeid = local_id;
    }

    // Store the xyz info of "layer" nodes
    fixed_layer_pt_xyz[ii].resize(3 * num_fixed_layer_node);
    PERIGEE_OMP_PARALLEL_FOR
    for(int nn=0; nn<num_fixed_layer_node; ++nn)
    {
      const int GID = fixed_layer_global_node[ii][nn];
      fixed_layer_pt_xyz[ii][3 * nn]     = vol_ctrlPts[3 * GID];
      fixed_layer_pt_xyz[ii][3 * nn + 1] = vol_ctrlPts[3 * GID + 1];
      fixed_layer_pt_xyz[ii][3 * nn + 2] = vol_ctrlPts[3 * GID + 2];
    }

    rotated_layer_global_node[ii] = rotated_layer_vien[ii];
    VEC_T::sort_unique_resize(rotated_layer_global_node[ii]);
    const int num_rotated_layer_node = VEC_T::get_size(rotated_layer_global_node[ii]);

    // Convert the GlobalNodeID in rotated_layer_vien to local indices in this ebc
    PERIGEE_OMP_PARALLEL_FOR
    for(int &nodeid : rotated_layer_vien[ii])
    {
      const int local_id = VEC_T::get_pos(rotated_layer_global_node[ii], nodeid);
      nodeid = local_id;
    }

    // Store the xyz info of "layer" nodes
    rotated_layer_pt_xyz[ii].resize(3 * num_rotated_layer_node);
    PERIGEE_OMP_PARALLEL_FOR
    for(int nn=0; nn<num_rotated_layer_node; ++nn)
    {
      const int GID = rotated_layer_global_node[ii][nn];
      rotated_layer_pt_xyz[ii][3 * nn]     = vol_ctrlPts[3 * GID];
      rotated_layer_pt_xyz[ii][3 * nn + 1] = vol_ctrlPts[3 * GID + 1];
      rotated_layer_pt_xyz[ii][3 * nn + 2] = vol_ctrlPts[3 * GID + 2];
    }
  }
}

ElemBC_3D_sliding_interface::~ElemBC_3D_sliding_interface()
{
  for(int ii=0; ii < num_interface_pair; ++ii)
  {
    fixed_face_id[ii].clear();
    fixed_part_tag[ii].clear();
    fixed_layer_vien[ii].clear();
    fixed_layer_global_node[ii].clear();
    fixed_layer_pt_xyz[ii].clear();
    rotated_face_id[ii].clear();
    rotated_layer_vien[ii].clear();
    rotated_layer_global_node[ii].clear();
    rotated_layer_pt_xyz[ii].clear();
  }

  fixed_face_id.clear();
  fixed_part_tag.clear();
  fixed_layer_vien.clear();
  fixed_layer_global_node.clear();
  fixed_layer_pt_xyz.clear();
  rotated_face_id.clear();
  rotated_layer_vien.clear();
  rotated_layer_global_node.clear();
  rotated_layer_pt_xyz.clear();
}

// EOF
