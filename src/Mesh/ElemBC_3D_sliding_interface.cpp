#include "ElemBC_3D_sliding_interface.hpp"

ElemBC_3D_sliding_interface::ElemBC_3D_sliding_interface( 
    const std::vector<std::string> &vtkfileList,
    const int &num_fixed_pt, const int &num_fixed_elem, const IIEN * const &VIEN, const int &elemtype )
: ElemBC_3D ( vtkfileList, elemtype )
{
  SYS_T::print_fatal_if(VEC_T::get_size(vtkfileList) != 2,
    "Error, ElemBC_3D_sliding_interface: There shold be 2 vtk files input, as the fixed/rotated interfaces.\n");
  
  fixed_face_id.resize(num_cell[0]);
  rotated_face_id.resize(num_cell[1]);

  if(elem_type == 501 || elem_type == 502)
    {
      TET_T::Tet4 * tetcell = new TET_T::Tet4();

      for(int ee{0}; ee < num_cell[0]; ++ee)
      {
        const int node_t[3] { get_ien(0, ee, 0), get_ien(0, ee, 1), get_ien(0, ee, 2) };

        const int node_t_gi[3] { get_global_node(0, node_t[0]),
                                 get_global_node(0, node_t[1]),
                                 get_global_node(0, node_t[2]) };

        const int cell_gi = get_global_cell(0, ee);

        const int tet_n[4] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3) };
        
        tetcell->reset(tet_n[0], tet_n[1], tet_n[2], tet_n[3]);

        fixed_face_id[ee] = tetcell->get_face_id(node_t_gi[0], node_t_gi[1], node_t_gi[2]);
      }

      for(int ee{0}; ee < num_cell[1]; ++ee)
      {
        const int node_t[3] { get_ien(1, ee, 0), get_ien(1, ee, 1), get_ien(1, ee, 2) };

        const int node_t_gi[3] { get_global_node(1, node_t[0]) + num_fixed_pt,
                                 get_global_node(1, node_t[1]) + num_fixed_pt,
                                 get_global_node(1, node_t[2]) + num_fixed_pt};

        global_cell[1][ee] += num_fixed_elem;
        const int cell_gi = get_global_cell(1, ee);

        const int tet_n[4] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3) };
        
        tetcell->reset(tet_n[0], tet_n[1], tet_n[2], tet_n[3]);

        rotated_face_id[ee] = tetcell->get_face_id(node_t_gi[0], node_t_gi[1], node_t_gi[2]);
      }
    }
    else if(elem_type == 601 || elem_type == 602)
    {
      HEX_T::Hex8 * hexcell = new HEX_T::Hex8();

      for(int ee{0}; ee < num_cell[0]; ++ee)
      {
        const int node_q[4] { get_ien(0, ee, 0), get_ien(0, ee, 1),
                              get_ien(0, ee, 2), get_ien(0, ee, 3) };
        
        const int node_q_gi[4] { get_global_node(0, node_q[0]),
                                 get_global_node(0, node_q[1]),
                                 get_global_node(0, node_q[2]),
                                 get_global_node(0, node_q[3]) };
        
        const int cell_gi = get_global_cell(0, ee);

        const int hex_n[8] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3),
                             VIEN->get_IEN(cell_gi, 4), VIEN->get_IEN(cell_gi, 5),
                             VIEN->get_IEN(cell_gi, 6), VIEN->get_IEN(cell_gi, 7) };
        
        hexcell->reset(hex_n[0], hex_n[1], hex_n[2], hex_n[3],
                       hex_n[4], hex_n[5], hex_n[6], hex_n[7]);

        fixed_face_id[ee] = hexcell->get_face_id(node_q_gi[0], node_q_gi[1], node_q_gi[2], node_q_gi[3]);
      }

      for(int ee{0}; ee < num_cell[1]; ++ee)
      {
        const int node_q[4] { get_ien(1, ee, 0), get_ien(1, ee, 1),
                              get_ien(1, ee, 2), get_ien(1, ee, 3) };
        
        const int node_q_gi[4] { get_global_node(1, node_q[0]) + num_fixed_pt,
                                 get_global_node(1, node_q[1]) + num_fixed_pt,
                                 get_global_node(1, node_q[2]) + num_fixed_pt,
                                 get_global_node(1, node_q[3]) + num_fixed_pt};
        
        global_cell[1][ee] += num_fixed_elem;
        const int cell_gi = get_global_cell(1, ee);

        const int hex_n[8] { VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3),
                             VIEN->get_IEN(cell_gi, 4), VIEN->get_IEN(cell_gi, 5),
                             VIEN->get_IEN(cell_gi, 6), VIEN->get_IEN(cell_gi, 7) };
        
        hexcell->reset(hex_n[0], hex_n[1], hex_n[2], hex_n[3],
                       hex_n[4], hex_n[5], hex_n[6], hex_n[7]);

        rotated_face_id[ee] = hexcell->get_face_id(node_q_gi[0], node_q_gi[1], node_q_gi[2], node_q_gi[3]);
      }
    }
    else
      SYS_T::print_fatal("Error: ElemBC_3D_sliding_interface, unknown element type.\n");
}

ElemBC_3D_sliding_interface::~ElemBC_3D_sliding_interface()
{
  fixed_face_id.clear();
  rotated_face_id.clear();
}

// EOF
