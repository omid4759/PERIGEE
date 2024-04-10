#include "ElemBC_3D_sliding_interface.hpp"

ElemBC_3D_sliding_interface::ElemBC_3D_sliding_interface( 
    const std::vector<std::string> &vtkfileList,
    const int &num_interface_pair_in,
    const int &num_fixed_pt, const int &num_fixed_elem,
    const std::vector<double> &vol_ctrlPts, const IIEN * const &VIEN, const int &elemtype )
: ElemBC_3D ( vtkfileList, elemtype ), num_interface_pair{num_interface_pair_in},
  fixed_face_id{std::vector<std::vector<int>> {}},  rotated_face_id{std::vector<std::vector<int>> {}},
  rotated_layer_vien{std::vector<std::vector<int>> {}}, rotated_layer_global_node{std::vector<std::vector<int>> {}},
  rotated_layer_pt_xyz{std::vector<std::vector<double>> {}}
{
  SYS_T::print_fatal_if(VEC_T::get_size(vtkfileList) != 2 * num_interface_pair,
    "Error, ElemBC_3D_sliding_interface: The number of interface file is wrong!\n");
  
  fixed_face_id.resize(num_interface_pair);
  rotated_face_id.resize(num_interface_pair);
  rotated_layer_vien.resize(num_interface_pair);
  for(int ii=0; ii<num_interface_pair; ++ii)
  {
    fixed_face_id[ii].resize(num_cell[ii]);
    rotated_face_id[ii].resize(num_cell[ii + num_interface_pair]);
    rotated_layer_vien[ii].resize(VIEN->get_nLocBas() * num_cell[ii + num_interface_pair]);
  }

  rotated_layer_global_node.resize(num_interface_pair);
  rotated_layer_pt_xyz.resize(num_interface_pair);

  if(elem_type == 501 || elem_type == 502)
  {
    TET_T::Tet4 * tetcell = new TET_T::Tet4();

    for(int ii=0; ii<num_interface_pair; ++ii)
    {
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
  }
  else if(elem_type == 601 || elem_type == 602)
  {
    HEX_T::Hex8 * hexcell = new HEX_T::Hex8();

    for(int ii=0; ii<num_interface_pair; ++ii)
    {
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
    rotated_face_id[ii].clear();
    rotated_layer_vien[ii].clear();
    rotated_layer_global_node[ii].clear();
    rotated_layer_pt_xyz[ii].clear();
  }

  fixed_face_id.clear();
  rotated_face_id.clear();
  rotated_layer_vien.clear();
  rotated_layer_global_node.clear();
  rotated_layer_pt_xyz.clear();
}

// EOF
