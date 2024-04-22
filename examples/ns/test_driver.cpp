// test program for sliding interface

#include "HDF5_Writer.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "ALocal_EBC_outflow.hpp"
#include "ALocal_WeakBC.hpp"
#include "ALocal_InflowBC.hpp"
#include "ALocal_Interface.hpp"
#include "QuadPts_UserDefined_Triangle.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "QuadPts_Gauss_Quad.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss_Hex.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Hex8.hpp"
#include "FEAElement_Hex27.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Quad4_3D_der0.hpp"
#include "FEAElement_Quad9_3D_der0.hpp"
#include "CVFlowRate_Unsteady.hpp"
#include "CVFlowRate_Linear2Steady.hpp"
#include "CVFlowRate_Cosine2Steady.hpp"
#include "GenBC_Resistance.hpp"
#include "GenBC_RCR.hpp"
#include "GenBC_Inductance.hpp"
#include "GenBC_Coronary.hpp"
#include "GenBC_Pressure.hpp"
#include "PLocAssem_VMS_NS_GenAlpha.hpp"
#include "PLocAssem_VMS_NS_GenAlpha_WeakBC.hpp"
#include "PGAssem_NS_FEM.hpp"
#include "PTime_NS_Solver.hpp"

int main(int argc, char *argv[])
{

  // Number of quadrature points for tets and triangles
  // Suggested values: 5 / 4 for linear, 17 / 13 for quadratic
  int nqp_tet = 5, nqp_tri = 4;
  
  // Estimate of the nonzero per row for the sparse matrix
  int nz_estimate = 300;

  double fluid_density = 1.065;
  double fluid_mu = 3.5e-2;
  double c_tauc = 1.0; // scaling factor for tau_c, take 0.0, 0.125, or 1.0
  double c_ct = 4.0; // C_T parameter for defining tau_M

  std::string lpn_file("lpn_pressure_input.txt");

  // back flow stabilization
  double bs_beta = 0.2;

  // generalized-alpha rho_inf
  double genA_rho_inf = 0.5;

  // part file location
  std::string part_file("part");


#if PETSC_VERSION_LT(3,19,0)
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
#else
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULLPTR);
#endif
  
  const PetscMPIInt rank = SYS_T::get_MPI_rank();

  SYS_T::print_perigee_art();

  // ===== Read Command Line Arguments =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionInt("-nz_estimate", nz_estimate);

  // ===== Record important solver options =====
  if(rank == 0)
  {
    hid_t cmd_file_id = H5Fcreate("solver_cmd.h5",
        H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

    cmdh5w->write_doubleScalar("fl_density", fluid_density);
    cmdh5w->write_doubleScalar("fl_mu", fluid_mu);
    cmdh5w->write_string("lpn_file", lpn_file);

    delete cmdh5w; H5Fclose(cmd_file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);

  FEANode * fNode = new FEANode(part_file, rank);

  // Local sub-domain's IEN array
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);

  // Global mesh info
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);

  // Mesh partition info
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);

  // Local sub-domain's element indices
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);

  // Local sub-domain's nodal bc
  ALocal_NBC * locnbc = new ALocal_NBC(part_file, rank);

  // Local sub-domain's inflow bc
  ALocal_InflowBC * locinfnbc = new ALocal_InflowBC(part_file, rank);

  // Local sub-domain's elemental bc
  ALocal_EBC * locebc = new ALocal_EBC_outflow(part_file, rank);

  // Local sub_domain's weak bc
  ALocal_WeakBC * locwbc = new ALocal_WeakBC(part_file, rank);
  locwbc -> print_info();

  ALocal_Interface* locitf = new ALocal_Interface(part_file, rank);
  locitf -> print_info();

  // Local sub-domain's nodal indices
  APart_Node * pNode = new APart_Node(part_file, rank);

  FEAElement * elementv = nullptr;
  FEAElement * elements = nullptr;
  FEAElement * elementvs = nullptr;
  IQuadPts * quadv = nullptr;
  IQuadPts * quads = nullptr;

  const int nqp_vol { nqp_tet };
  const int nqp_sur { nqp_tri };

  elementv = new FEAElement_Tet4( nqp_vol ); // elem type 501
  elements = new FEAElement_Triangle3_3D_der0( nqp_sur );
  elementvs = new FEAElement_Tet4( nqp_sur );
  quadv = new QuadPts_Gauss_Tet( nqp_vol );
  quads = new QuadPts_Gauss_Triangle( nqp_sur );

  // ===== Generate a sparse matrix for the enforcement of essential BCs
  Matrix_PETSc * pmat = new Matrix_PETSc(pNode, locnbc);

  pmat->gen_perm_bc(pNode, locnbc);

  // ===== Generalized-alpha =====
  SYS_T::commPrint("===> Setup the Generalized-alpha time scheme.\n");

  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(
      genA_rho_inf, false );

  IPLocAssem * locAssem_ptr = nullptr;
  locAssem_ptr = new PLocAssem_VMS_NS_GenAlpha(
      tm_galpha_ptr, elementv->get_nLocBas(),
      quadv->get_num_quadPts(), elements->get_nLocBas(),
      fluid_density, fluid_mu, bs_beta, GMIptr->get_elemType(), c_ct, c_tauc );

  IGenBC * gbc = nullptr;
  gbc = new GenBC_Pressure( lpn_file, 0.0 );

  PGAssem_NS_FEM * gloAssem_ptr = new PGAssem_NS_FEM( locAssem_ptr, elements, quads,
      GMIptr, locElem, locIEN, pNode, locnbc, locebc, gbc, nz_estimate );

  PDNSolution * sol = new PDNSolution_NS( pNode, 0 );

  IQuadPts * free_quad = new QuadPts_UserDefined_Triangle(nqp_sur);

  gloAssem_ptr->Interface_G(0, 0, sol, locAssem_ptr, elementvs, elementv, elements, quads, free_quad, locIEN, fNode, locitf);

  delete fNode; delete locIEN; delete GMIptr; delete PartBasic;
  delete locElem; delete locnbc; delete locebc; delete locwbc; delete pNode; delete locinfnbc; delete locitf;
  delete tm_galpha_ptr; delete pmat; delete elementv; delete elements; delete elementvs;
  delete quads; delete quadv;

  delete locAssem_ptr; delete gbc; delete gloAssem_ptr;
  delete sol; delete free_quad;

  SYS_T::commPrint("Successfully end.\n");

  PetscFinalize();
  return EXIT_SUCCESS;
}