#include "PDNSolution_Wall_Disp.hpp"

PDNSolution_Wall_Disp::PDNSolution_Wall_Disp( 
    const APart_Node * const &pNode,
    const int &type, const bool &isprint )
: PDNSolution( pNode, 3 ), is_print( isprint )
{
  switch( type )
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: in PDNSolution_Wall_Disp, No such type of initial condition. \n");
      break;
  }
}


// Initialize from vtu file with displacement fields
PDNSolution_Wall_Disp::PDNSolution_Wall_Disp(
    const APart_Node * const &pNode,
    const std::string &init_vtu,
    const std::string &disp_name,
    const bool &isprint )
: PDNSolution( pNode, 3 ), is_print( isprint )
{
  SYS_T::file_check( init_vtu );

  std::vector<double> init_disp = TET_T::read_double_PointData( init_vtu, disp_name, 3 );

  double value[3];
  const int nlocalnode = pNode->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode->get_node_loc(ii) * 3;
    const int location[3] = { pos, pos + 1, pos +2 };

    // Access original node idx prior to old2new mapping
    const int original_idx = pNode->get_node_loc_original(ii);

    value[0] = init_disp[3 * original_idx + 0]; 
    value[1] = init_disp[3 * original_idx + 1]; 
    value[2] = init_disp[3 * original_idx + 2]; 

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);

  GhostUpdate();

  if(is_print)
  {
    SYS_T::commPrint("===> Initial solution: Read the following fields from %s\n", init_vtu.c_str());
    SYS_T::commPrint("                       disp_x = %s \n", disp_name.c_str());
    SYS_T::commPrint("                       disp_y = %s \n", disp_name.c_str());
    SYS_T::commPrint("                       disp_z = %s \n", disp_name.c_str());
  }
}


PDNSolution_Wall_Disp::~PDNSolution_Wall_Disp()
{}


void PDNSolution_Wall_Disp::Init_zero( const APart_Node * const &pNode_ptr )
{
  const double value[3] = {0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode_ptr->get_node_loc(ii) * 3;
    const int location[3] = { pos, pos + 1, pos + 2 };

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);

  GhostUpdate();

  if(is_print)
  {
    SYS_T::commPrint("===> Initial solution: disp_x = 0.0 \n");
    SYS_T::commPrint("                       disp_y = 0.0 \n");
    SYS_T::commPrint("                       disp_z = 0.0 \n");
  }
}

// EOF
