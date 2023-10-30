#include "PDNSolution_V.hpp"

PDNSolution_V::PDNSolution_V( const APart_Node * const &pNode,
    const int &type, const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode ), sol_name( in_name ), is_print( isprint )
{
  SYS_T::print_fatal_if( pNode->get_dof() != 3, "Error: PDNSolution_V : the APart_Node gives wrong dof number. \n");
  
  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_V: No such type of initial condition. \n");
      break;
  }
}

PDNSolution_V::PDNSolution_V( const APart_Node * const &pNode,
    const FEANode * const &fNode,
    const ALocal_InflowBC * const &infbc,
    const int &type, const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode ), sol_name( in_name ), is_print( isprint )
{
  SYS_T::print_fatal_if( pNode->get_dof() != 3, "Error: PDNSolution_V : the APart_Node gives wrong dof number. \n");
  
  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    case 1:
      Init_flow_parabolic( pNode, fNode, infbc );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_V: No such type of initial condition. \n");
      break;
  }
}

PDNSolution_V::PDNSolution_V( const APart_Node * const &pNode,
    const FEANode * const &fNode,
    const double &rho,
    const double &vis_mu,
    const int &type, const bool &isprint,
    const std::string &in_name )
: PDNSolution( pNode ), sol_name( in_name ), is_print( isprint )
{
  SYS_T::print_fatal_if( pNode->get_dof() != 3, "Error: PDNSolution_V : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    case 2:
      Init_flow_womersley( pNode, fNode, rho, vis_mu );
      break;
    case 3:
      Init_flow_womersley_dot( pNode, fNode, rho, vis_mu );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_V: No such type of initial condition. \n");
      break;
  }
}

void PDNSolution_V::Init_zero( const APart_Node * const &pNode )
{
  double value[3] = {0.0, 0.0, 0.0};

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode -> get_node_loc(ii) * 3;
    const int location[3] = { pos, pos + 1, pos + 2 };

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
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

void PDNSolution_V::Init_flow_parabolic( const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const ALocal_InflowBC * const &infbc )
{
  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode_ptr->get_node_loc(ii) * 3;
    const int location[3] = { pos, pos + 1, pos + 2 };
    const double value[3] = {0.0, 0.0, 0.0};

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  const int num_nbc = infbc -> get_num_nbc();

  for(int nbc_id = 0; nbc_id < num_nbc; ++nbc_id)
  {
    const double vmax = 2.0 / infbc->get_fularea( nbc_id );
    const double out_nx = infbc->get_outvec( nbc_id ).x();
    const double out_ny = infbc->get_outvec( nbc_id ).y();
    const double out_nz = infbc->get_outvec( nbc_id ).z();

    // If this sub-domain contains inflow nodes, set their values based on the
    // parabolic flow profile
    if( infbc->get_Num_LD( nbc_id ) > 0)
    {
      for(int ii=0; ii<nlocalnode; ++ii)
      {
        if( infbc->is_inLDN( nbc_id, pNode_ptr->get_node_loc(ii) ) )
        {
          const int pos = pNode_ptr->get_node_loc(ii) * 3;
          const int location[3] = { pos, pos + 1, pos + 2 };

          const Vector_3 pt = fNode_ptr -> get_ctrlPts_xyz(ii);
          const double r =  infbc -> get_radius( nbc_id, pt );
          const double vel = vmax * (1.0 - r*r);

          const double value[3] = { vel * out_nx, vel * out_ny, vel * out_nz };

          VecSetValues(solution, 3, location, value, INSERT_VALUES);
        }
      }
    }
  }

  Assembly_GhostUpdate();

  if( is_print )
  {
    std::ostringstream ss;
    ss<<"===> Initial "<<sol_name<<" solution vector: \n";
    SYS_T::commPrint(ss.str().c_str());
    for(int nbc_id=0; nbc_id < num_nbc; ++nbc_id)
    {
      SYS_T::commPrint("     -- nbc_id = %d \n", nbc_id);
      SYS_T::commPrint("        max speed %e.\n", 2.0 / infbc->get_fularea( nbc_id ) );
      SYS_T::commPrint("        active area is %e.\n", infbc->get_actarea(nbc_id) );
      SYS_T::commPrint("        full area is %e.\n", infbc->get_fularea(nbc_id) );
      SYS_T::commPrint("        outward normal direction [%e %e %e].\n",
          infbc->get_outvec( nbc_id ).x(), infbc->get_outvec( nbc_id ).y(), infbc->get_outvec( nbc_id ).z() );
    }
  }
}

void PDNSolution_V::Init_flow_womersley(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const double &rho,
    const double &vis_mu )
{
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode_ptr->get_node_loc(ii) * 3;
    const int location[3] = { pos, pos + 1, pos + 2 };
    const double value[3] = {0.0, 0.0, 0.0};

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  const double R     = 0.3;                                                  // pipe radius
  const double omega = MATH_T::PI * 2.0 / 1.1;                               // freqency
  const std::complex<double> i1(0.0, 1.0);
  const std::complex<double> i1_1d5(-0.707106781186547, 0.707106781186547);
  const auto Omega   = std::sqrt(rho * omega / vis_mu) * R;                  // womersley number
  const auto Lambda  = i1_1d5 * Omega;

  const double k0 = -21.0469;                                                // mean pressure gradient
  const std::complex<double> B1(-4.926286624202966e3, -4.092542965905093e3); // pressure Fourier coeff
  const std::complex<double> c1(8.863128942479001e2,   2.978553160539686e1); // wave speed
  const std::complex<double> G1(0.829733473284180,      -0.374935589823809); // elasticity factor

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode_ptr->get_node_loc(ii) * 3;
    const int location[3] = { pos, pos + 1, pos + 2 };

    const double x  = fNode_ptr->get_ctrlPts_x(ii);
    const double y  = fNode_ptr->get_ctrlPts_y(ii);
    const double z  = fNode_ptr->get_ctrlPts_z(ii);
    const double r  = std::sqrt(x*x + y*y);
    const auto   xi = Lambda * r / R;

    const auto bes0_xi     = sp_bessel::besselJ(0, xi);
    const auto bes1_xi     = sp_bessel::besselJ(1, xi);
    const auto bes0_Lambda = sp_bessel::besselJ(0, Lambda);

    // pressure
    // const double pres = k0 * z + std::real( B1 * exp(-i1*omega*z/c1) );

    // axial velocity
    const double w = k0 * (x*x + y*y - R*R) / (4.0*vis_mu)
        + std::real( B1 / (rho * c1) * (1.0 - G1 * bes0_xi / bes0_Lambda) * exp(-i1*omega*z/c1) );

    // radial velo
    const double vr = std::real( i1 * omega * R * B1 / ( 2.0 * rho * c1 * c1 )
        * ( r / R - 2.0 * G1 * bes1_xi / (Lambda * bes0_Lambda) ) * exp(-i1*omega*z/c1) );

    // polar to cartesian transformation
    const double theta = std::atan2(y, x);

    const double value[3] = { vr * std::cos(theta), vr * std::sin(theta), w };

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();
}

void PDNSolution_V::Init_flow_womersley_dot(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const double &rho,
    const double &vis_mu )
{
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode_ptr->get_node_loc(ii) * 3;
    const int location[3] = { pos, pos + 1, pos + 2 };
    const double value[3] = {0.0, 0.0, 0.0};

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  const double R     = 0.3;                                                  // pipe radius
  const double omega = MATH_T::PI * 2.0 / 1.1;                               // freqency
  const std::complex<double> i1(0.0, 1.0);
  const std::complex<double> i1_1d5(-0.707106781186547, 0.707106781186547);
  const auto Omega   = std::sqrt(rho * omega / vis_mu) * R;                  // womersley number
  const auto Lambda  = i1_1d5 * Omega;

  // const double k0 = -21.0469;                                                // mean pressure gradient
  const std::complex<double> B1(-4.926286624202966e3, -4.092542965905093e3); // pressure Fourier coeff
  const std::complex<double> c1(8.863128942479001e2,   2.978553160539686e1); // wave speed
  const std::complex<double> G1(0.829733473284180,      -0.374935589823809); // elasticity factor

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    const int pos = pNode_ptr->get_node_loc(ii) * 3;
    const int location[3] = { pos, pos + 1, pos + 2 };

    const double x  = fNode_ptr->get_ctrlPts_x(ii);
    const double y  = fNode_ptr->get_ctrlPts_y(ii);
    const double z  = fNode_ptr->get_ctrlPts_z(ii);
    const double r  = std::sqrt(x*x + y*y);
    const auto   xi = Lambda * r / R;

    const auto bes0_xi     = sp_bessel::besselJ(0, xi);
    const auto bes1_xi     = sp_bessel::besselJ(1, xi);
    const auto bes0_Lambda = sp_bessel::besselJ(0, Lambda);

    // dot pressure
    // const double dot_pres = std::real( i1 * omega * B1 * exp(-i1*omega*z/c1) );

    // dot axial velocity
    const double dot_w = std::real( i1 * omega * B1 / (rho * c1)
        * (1.0 - G1 * bes0_xi / bes0_Lambda) * exp(-i1*omega*z/c1) );

    // dot radial velocity
    const double dot_vr = std::real( -omega * omega * R * B1 / ( 2.0 * rho * c1 * c1 )
        * ( r / R - 2.0 * G1 * bes1_xi / (Lambda * bes0_Lambda) ) * exp(-i1*omega*z/c1) );

    // polar to cartesian transformation
    const double theta = std::atan2(y, x);

    const double value[3] = { dot_vr * std::cos(theta), dot_vr * std::sin(theta), dot_w };

    VecSetValues(solution, 3, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();
}

// EOF
