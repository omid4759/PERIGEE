#include "PDNSolution_P.hpp"

PDNSolution_P::PDNSolution_P( const APart_Node * const &pNode,
    const int &type, const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode ), sol_name( in_name ), is_print( isprint )
{
  SYS_T::print_fatal_if( pNode->get_dof() != 1, "Error: PDNSolution_P : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_P: No such type of initial condition. \n");
      break;
  }
}

PDNSolution_P::PDNSolution_P( const APart_Node * const &pNode,
    const FEANode * const &fNode,
    const int &type, const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode ), sol_name( in_name ), is_print( isprint )
{
  SYS_T::print_fatal_if( pNode->get_dof() != 1, "Error: PDNSolution_P : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    case 2:
      Init_pres_womersley( pNode, fNode );
      break;
    case 3:
      Init_pres_womersley_dot( pNode, fNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_P: No such type of initial condition. \n");
      break;
  }
}

void PDNSolution_P::Init_zero( const APart_Node * const &pNode )
{
  const double val = 0.0;

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode -> get_node_loc(ii);
    VecSetValue(solution, pos, val, INSERT_VALUES);
  }
  
  Assembly_GhostUpdate();

  if( is_print )
  {
    std::ostringstream ss;
    ss<<"===> Initial "<<sol_name<<" solution vector: \n";
    SYS_T::commPrint(ss.str().c_str());
    SYS_T::commPrint("     val_x = 0.0 \n");
    SYS_T::commPrint("     val_y = 0.0 \n");
    SYS_T::commPrint("     val_z = 0.0 \n");
  }
}

void PDNSolution_P::Init_pres_womersley( const APart_Node * const &pNode,
    const FEANode * const &fNode )
{
  const int nlocalnode = pNode->get_nlocalnode();
  
  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode -> get_node_loc(ii);
    const double val = 0.0;
    VecSetValue(solution, pos, val, INSERT_VALUES);
  }

  const double omega = MATH_T::PI * 2.0 / 1.1;                               // freqency
  const std::complex<double> i1(0.0, 1.0);

  const double k0 = -21.0469;                                                // mean pressure gradient
  const std::complex<double> B1(-4.926286624202966e3, -4.092542965905093e3); // pressure Fourier coeff
  const std::complex<double> c1(8.863128942479001e2,   2.978553160539686e1); // wave speed

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode->get_node_loc(ii);
    const double z  = fNode->get_ctrlPts_z(ii);

    // pressure
    const double pres = k0 * z + std::real( B1 * exp(-i1*omega*z/c1) );

    const double val = pres;
    VecSetValue(solution, pos, val, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();
}

void PDNSolution_P::Init_pres_womersley_dot( const APart_Node * const &pNode,
    const FEANode * const &fNode )
{
  const int nlocalnode = pNode->get_nlocalnode();

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode -> get_node_loc(ii);
    const double val = 0.0;
    VecSetValue(solution, pos, val, INSERT_VALUES);
  }
                           
  const double omega = MATH_T::PI * 2.0 / 1.1;                               // freqency
  const std::complex<double> i1(0.0, 1.0);

  // const double k0 = -21.0469;                                                // mean pressure gradient
  const std::complex<double> B1(-4.926286624202966e3, -4.092542965905093e3); // pressure Fourier coeff
  const std::complex<double> c1(8.863128942479001e2,   2.978553160539686e1); // wave speed

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode->get_node_loc(ii);
    const double z  = fNode->get_ctrlPts_z(ii);

    // dot pressure
    const double dot_pres = std::real( i1 * omega * B1 * exp(-i1*omega*z/c1) );

    const double val = dot_pres;
    VecSetValue(solution, pos, val, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();
}


// EOF
