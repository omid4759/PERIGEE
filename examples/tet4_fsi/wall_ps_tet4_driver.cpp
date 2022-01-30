// ============================================================================
// wall_ps_tet4_driver.cpp
//
// Wall mechanics solver for generating the prestress field.
//
// Date: Jan 28 2022
// ============================================================================
#include "HDF5_Tools.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node_FSI.hpp"
#include "ALocal_Elem.hpp"


#include "ALocal_EBC_outflow.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Tet4.hpp"
#include "CVFlowRate_Unsteady.hpp"
#include "CVFlowRate_Linear2Steady.hpp"
#include "CVFlowRate_Steady.hpp"
#include "GenBC_Resistance.hpp"
#include "GenBC_RCR.hpp"
#include "GenBC_Inductance.hpp"
#include "GenBC_Coronary.hpp"
#include "GenBC_Pressure.hpp"
#include "MaterialModel_NeoHookean_M94_Mixed.hpp"
#include "MaterialModel_NeoHookean_Incompressible_Mixed.hpp"
#include "PLocAssem_2x2Block_Tet4_ALE_VMS_NS_GenAlpha.hpp"
#include "PLocAssem_2x2Block_Tet4_VMS_Incompressible.hpp"
#include "PLocAssem_2x2Block_Tet4_VMS_Hyperelasticity.hpp"
#include "PLocAssem_Tet4_FSI_Mesh_Elastostatic.hpp"
#include "PLocAssem_Tet4_FSI_Mesh_Laplacian.hpp"
#include "PGAssem_FSI.hpp"
#include "PGAssem_Mesh.hpp"
#include "PTime_FSI_Solver.hpp"

#include "PDNTimeStep.hpp"
#include "PETSc_Tools.hpp"

int main( int argc, char *argv[] )
{
  // solution file name to be loaded for prestressing
  std::string restart_velo_name = "SOL_velo_re";
  std::string restart_pres_name = "SOL_pres_re";

  // (Pseudo-) time integration parameters
  double genA_rho_inf = 0.0;
  bool is_backward_Euler = true;
  const bool is_load_ps = false;

  // Estimate of num nonzeros per row for the sparse tangent matrix
  int nz_estimate = 300;

  // Prestress tolerance
  double prestress_disp_tol = 1.0e-6;

  // Nonlinear solver parameters
  double nl_rtol = 1.0e-3;           // convergence criterion relative tolerance
  double nl_atol = 1.0e-6;           // convergence criterion absolute tolerance
  double nl_dtol = 1.0e3;            // divergence criterion
  int    nl_maxits = 20;             // maximum number if nonlinear iterations
  int    nl_refreq = 4;              // frequency of tangent matrix renewal
  int    nl_threshold = 4;           // threshold of tangent matrix renewal

  // Time stepping parameters
  double initial_time = 0.0;         // time of initial condition
  double initial_step = 0.1;         // time step size
  int    initial_index = 0;          // index of initial condition
  double final_time = 1.0;           // end time of simulation
  bool   is_record_sol = false;      // bool flag to decide if one wants to record the solution
  std::string sol_bName("PS_");      // base name of the solution file
  int    ttan_renew_freq = 1;        // frequency of tangent matrix renewal
  int    sol_record_freq = 1;        // frequency for recording the solution

  // We assume that a 3D solver has been called (to generate the wall traction)
  // and a suite of command line arguments has been saved to disk
  hid_t solver_cmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( solver_cmd_file );

  const int nqp_tet       = cmd_h5r -> read_intScalar(    "/", "nqp_tet");
  const int nqp_tri       = cmd_h5r -> read_intScalar(    "/", "nqp_tri");
  const double sl_nu      = cmd_h5r -> read_doubleScalar( "/", "sl_nu");

  delete cmd_h5r; H5Fclose(solver_cmd_file);

  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5_Reader * pcmd_h5r = new HDF5_Reader( prepcmd_file );

  const std::string part_v_file = pcmd_h5r -> read_string(    "/", "part_file_v" );
  const std::string part_p_file = pcmd_h5r -> read_string(    "/", "part_file_p" );
  const int fsiBC_type        = pcmd_h5r -> read_intScalar( "/", "fsiBC_type" );

  delete pcmd_h5r; H5Fclose(prepcmd_file);

  // Initialize PETSc
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::print_fatal_if( fsiBC_type != 2, "Error: fsiBC_type should be 2.\n" );

  // Clean potentially pre-existing hdf5 files of prestress saved in the folder
  // named as prestress
  if(rank == 0 )
  {
    if( SYS_T::directory_exist("prestress") )
    {
      std::cout<<"Clean the folder prestress.\n";
      SYS_T::execute("rm -rf prestress");
    }

    SYS_T::execute("mkdir prestress");
  }

  SYS_T::GetOptionString("-restart_velo_name",   restart_velo_name);
  SYS_T::GetOptionString("-restart_pres_name",   restart_pres_name);
  SYS_T::GetOptionReal(  "-rho_inf",             genA_rho_inf);
  SYS_T::GetOptionBool(  "-is_backward_Euler",   is_backward_Euler);
  SYS_T::GetOptionInt(   "-nz_estimate",         nz_estimate);
  SYS_T::GetOptionReal(  "-prestress_disp_tol",  prestress_disp_tol);
  SYS_T::GetOptionReal(  "-nl_rtol",             nl_rtol);
  SYS_T::GetOptionReal(  "-nl_atol",             nl_atol);
  SYS_T::GetOptionReal(  "-nl_dtol",             nl_dtol);
  SYS_T::GetOptionInt(   "-nl_maxits",           nl_maxits);
  SYS_T::GetOptionInt(   "-nl_refreq",           nl_refreq);
  SYS_T::GetOptionInt(   "-nl_threshold",        nl_threshold);
  SYS_T::GetOptionBool(  "-is_backward_Euler",   is_backward_Euler);
  SYS_T::GetOptionReal(  "-init_time",           initial_time);
  SYS_T::GetOptionReal(  "-fina_time",           final_time);
  SYS_T::GetOptionReal(  "-init_step",           initial_step);
  SYS_T::GetOptionInt(   "-init_index",          initial_index);
  SYS_T::GetOptionInt(   "-ttan_freq",           ttan_renew_freq);
  SYS_T::GetOptionBool(  "-is_record_sol",       is_record_sol);
  SYS_T::GetOptionInt(   "-sol_rec_freq",        sol_record_freq);

  // ===== Print Command Line Arguments =====
  SYS_T::cmdPrint(      "part_v_file:",          part_v_file);
  SYS_T::cmdPrint(      "part_p_file:",          part_p_file);
  SYS_T::cmdPrint(       "-prestress_disp_tol:", prestress_disp_tol);
  SYS_T::cmdPrint(       "-nl_rtol:",            nl_rtol);
  SYS_T::cmdPrint(       "-nl_atol:",            nl_atol);
  SYS_T::cmdPrint(       "-nl_dtol:",            nl_dtol);
  SYS_T::cmdPrint(       "-nl_maxits:",          nl_maxits);
  SYS_T::cmdPrint(       "-nl_refreq:",          nl_refreq);
  SYS_T::cmdPrint(       "-nl_threshold:",       nl_threshold);

  if( is_backward_Euler )
    SYS_T::commPrint(    "-is_backward_Euler: true \n");
  else
    SYS_T::cmdPrint(     "-rho_inf:",            genA_rho_inf);

  SYS_T::cmdPrint(       "-init_time:",          initial_time);
  SYS_T::cmdPrint(       "-init_step:",          initial_step);
  SYS_T::cmdPrint(       "-init_index:",         initial_index);
  SYS_T::cmdPrint(       "-fina_time:",          final_time);
  SYS_T::cmdPrint(       "-ttan_freq:",          ttan_renew_freq);

  if( is_record_sol )
    SYS_T::cmdPrint(     "-sol_rec_freq:",       sol_record_freq);
  else
    SYS_T::commPrint(    "-is_record_sol: false \n");

  // ====== Data for Analysis ======
  FEANode * fNode = new FEANode(part_v_file, rank);

  ALocal_IEN * locIEN_v = new ALocal_IEN(part_v_file, rank);

  ALocal_IEN * locIEN_p = new ALocal_IEN(part_p_file, rank);

  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_v_file, rank);

  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_v_file, rank);

  ALocal_Elem * locElem = new ALocal_Elem(part_v_file, rank);

  APart_Node * pNode_v = new APart_Node_FSI(part_v_file, rank);

  APart_Node * pNode_p = new APart_Node_FSI(part_p_file, rank);

  ALocal_EBC * locebc_v = new ALocal_EBC_outflow(part_v_file, rank);

  ALocal_EBC * locebc_p = new ALocal_EBC( part_p_file, rank );

  ALocal_EBC * mesh_locebc = new ALocal_EBC(part_v_file, rank, "/mesh_ebc");  

  ALocal_NodalBC * locnbc_v = new ALocal_NodalBC(part_v_file, rank, "/nbc/MF");

  ALocal_NodalBC * locnbc_p = new ALocal_NodalBC(part_p_file, rank, "/nbc/MF");

  ALocal_NodalBC * mesh_locnbc = new ALocal_NodalBC(part_v_file, rank, "/mesh_nbc/MF");

  Prestress_solid * ps_data = new Prestress_solid(locElem, nqp_tet, rank, is_load_ps, "prestress");  

  SYS_T::commPrint("===> Mesh HDF5 files are read from disk.\n");

  // Group APart_Node and ALocal_NodalBC into a vector
  std::vector<APart_Node *> pNode_list { pNode_v, pNode_p };

  std::vector<ALocal_NodalBC *> locnbc_list { locnbc_v, locnbc_p };

  std::vector<APart_Node *> pNode_m_list { pNode_v };

  std::vector<ALocal_NodalBC *> locnbc_m_list { mesh_locnbc };

  // ===== Basic Checking =====
  SYS_T::print_fatal_if( size!= PartBasic->get_cpu_size(),
      "Error: Assigned CPU number does not match the partition. \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for FEM analysis. \n", size);

  // ===== Quadrature rules and FEM container =====
  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );

  SYS_T::commPrint("===> Setup element container. \n");
  FEAElement * elementv = new FEAElement_Tet4( nqp_tet );
  FEAElement * elements = new FEAElement_Triangle3_3D_der0( nqp_tri );

  // ===== Generate the IS for pres and velo =====
  const int idx_v_start = HDF5_T::read_intScalar( SYS_T::gen_partfile_name(part_v_file, rank).c_str(), "/DOF_mapper", "start_idx" );
  const int idx_p_start = HDF5_T::read_intScalar( SYS_T::gen_partfile_name(part_p_file, rank).c_str(), "/DOF_mapper", "start_idx" );

  const int idx_v_len = pNode_v->get_dof() * pNode_v -> get_nlocalnode();
  const int idx_p_len = pNode_p->get_dof() * pNode_p -> get_nlocalnode();

  PetscInt * is_array_velo = new PetscInt[ idx_v_len ];
  for(int ii=0; ii<idx_v_len; ++ii) is_array_velo[ii] = idx_v_start + ii;

  PetscInt * is_array_pres = new PetscInt[ idx_p_len ];
  for(int ii=0; ii<idx_p_len; ++ii) is_array_pres[ii] = idx_p_start + ii;

  IS is_velo, is_pres;
  ISCreateGeneral(PETSC_COMM_WORLD, idx_v_len, is_array_velo, PETSC_COPY_VALUES, &is_velo);
  ISCreateGeneral(PETSC_COMM_WORLD, idx_p_len, is_array_pres, PETSC_COPY_VALUES, &is_pres);

  delete [] is_array_velo; is_array_velo = nullptr;
  delete [] is_array_pres; is_array_pres = nullptr;
  // ================================================================

  // ===== Generate a sparse matrix for strong enforcement of essential BC
  Matrix_PETSc * pmat = new Matrix_PETSc( idx_v_len + idx_p_len );
  pmat -> gen_perm_bc( pNode_list, locnbc_list );

  Matrix_PETSc * mmat = new Matrix_PETSc( pNode_v -> get_nlocalnode() * pNode_v -> get_dof() );
  mmat -> gen_perm_bc( pNode_m_list, locnbc_m_list );

  // ===== Generate the generalized-alpha method
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");

  TimeMethod_GenAlpha * tm_galpha_ptr = nullptr;

  if( is_backward_Euler )
    tm_galpha_ptr = new TimeMethod_GenAlpha( 1.0, 1.0, 1.0 );
  else
    tm_galpha_ptr = new TimeMethod_GenAlpha( genA_rho_inf, false );

  tm_galpha_ptr->print_info();

  // ===== Local assembly =====
  IMaterialModel * matmodel       = nullptr;
  IPLocAssem_2x2Block * locAssem_solid_ptr = nullptr;

  if( sl_nu == 0.5 )
  {
    matmodel = new MaterialModel_NeoHookean_Incompressible_Mixed( "material_model.h5" );

    locAssem_solid_ptr = new PLocAssem_2x2Block_Tet4_VMS_Incompressible(
        matmodel, tm_galpha_ptr, elementv -> get_nLocBas(), elements->get_nLocBas() );
  }
  else
  {
    matmodel = new MaterialModel_NeoHookean_M94_Mixed( "material_model.h5" );

    locAssem_solid_ptr = new PLocAssem_2x2Block_Tet4_VMS_Hyperelasticity(
        matmodel, tm_galpha_ptr, elementv -> get_nLocBas(), elements->get_nLocBas() );
  }

  // ===== Initial conditions =====
  PDNSolution * velo = new PDNSolution_V(pNode_v, 0, true, "velo");
  PDNSolution * disp = new PDNSolution_V(pNode_v, 0, true, "disp");
  PDNSolution * pres = new PDNSolution_P(pNode_p, 0, true, "pres");

  PDNSolution * dot_velo = new PDNSolution_V(pNode_v, 0, true, "dot_velo");
  PDNSolution * dot_disp = new PDNSolution_V(pNode_v, 0, true, "dot_disp");
  PDNSolution * dot_pres = new PDNSolution_P(pNode_p, 0, true, "dot_pres");

  // Read sol file
  SYS_T::file_check(restart_velo_name.c_str());
  velo->ReadBinary(restart_velo_name.c_str());

  SYS_T::file_check(restart_pres_name.c_str());
  pres->ReadBinary(restart_pres_name.c_str());

  // Read dot_sol file
  std::string restart_dot_velo_name = "dot_";
  restart_dot_velo_name.append(restart_velo_name);
  SYS_T::file_check(restart_dot_velo_name.c_str());
  dot_velo->ReadBinary(restart_dot_velo_name.c_str());

  std::string restart_dot_pres_name = "dot_";
  restart_dot_pres_name.append(restart_pres_name);
  SYS_T::file_check(restart_dot_pres_name.c_str());
  dot_pres->ReadBinary(restart_dot_pres_name.c_str());

  SYS_T::commPrint("===> Read sol from disk as a restart run: \n");
  SYS_T::commPrint("     restart_velo_name:     %s \n", restart_velo_name.c_str());
  SYS_T::commPrint("     restart_dot_velo_name: %s \n", restart_dot_velo_name.c_str());
  SYS_T::commPrint("     restart_pres_name:     %s \n", restart_pres_name.c_str());
  SYS_T::commPrint("     restart_dot_pres_name: %s \n", restart_dot_pres_name.c_str());

  // ===== Time step info =====
  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

























  PetscFinalize();
  return EXIT_SUCCESS;
}


// EOF 