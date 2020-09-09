// ==================================================================
// vis_ns.cpp
//
// Visualization driver for the NS equation solver.
// 
// Author: Ju Liu, liujuy@gmail.com
// Date Created: Feb. 12 2020
// ==================================================================
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "QuadPts_vis_tet10_v2.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "VisDataPrep_NS.hpp"
#include "VTK_Writer_NS.hpp"

// ==== WOMERSLEY CHANGES BEGIN ====
#include "PostVectSolution.hpp"
#include "Post_error_ns.hpp"
// ==== WOMERSLEY CHANGES END ====

int main( int argc, char * argv[] )
{
  const std::string element_part_file = "epart.h5";
  const std::string anode_mapping_file = "node_mapping.h5";
  const std::string pnode_mapping_file = "post_node_mapping.h5";
  const std::string part_file="postpart";
  const int dof = 4;
  
  std::string sol_bname("SOL_");
  std::string out_bname = sol_bname;
  int time_start = 0;
  int time_step = 1;
  int time_end = 1;
  bool isXML = true;
  bool isRestart = false;

  // Read analysis code parameter if the solver_cmd.h5 exists
  hid_t prepcmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );

  double dt = cmd_h5r -> read_doubleScalar("/","init_step");

  delete cmd_h5r; H5Fclose(prepcmd_file);

  // ===== Initialize the MPI run =====
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);

  PetscMPIInt rank = SYS_T::get_MPI_rank();
  PetscMPIInt size = SYS_T::get_MPI_size();

  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionReal("-dt", dt);
  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionBool("-restart", isRestart);

  // ==== WOMERSLEY CHANGES BEGIN ==== 
  int manu_sol_time = 0;
  SYS_T::GetOptionInt("-manu_sol_time", manu_sol_time);
  // ==== WOMERSLEY CHANGES END ====
  
  SYS_T::commPrint("=== Command line arguments ===\n");
  SYS_T::cmdPrint("-sol_bname:", sol_bname);
  SYS_T::cmdPrint("-out_bname:", out_bname);
  SYS_T::cmdPrint("-time_start:", time_start);
  SYS_T::cmdPrint("-time_step:", time_step);
  SYS_T::cmdPrint("-time_end:", time_end);
  SYS_T::cmdPrint("-dt:",dt);
  // ==== WOMERSLEY CHANGES BEGIN ==== 
  SYS_T::cmdPrint("-manu_sol_time:", manu_sol_time);
  // ==== WOMERSLEY CHANGES END ====

  if(isXML) SYS_T::commPrint("-xml: true \n");
  else SYS_T::commPrint("-xml: false \n");

  if(isRestart) SYS_T::commPrint("-restart: true \n");
  else SYS_T::commPrint("-restart: false \n");
  SYS_T::commPrint("==============================\n");
  
  // Clean the visualization files if not restart
  if( !isRestart )
  {
    int sysret = system("rm -rf *_p*.vtu");
    SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
    sysret = system("rm -rf *.pvtu");
    SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
    sysret = system("rm -rf *_.pvd");
    SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  }
  
  FEANode * fNode = new FEANode(part_file, rank);
  
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);
  
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);
  
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
 
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);
  
  APart_Node * pNode = new APart_Node(part_file, rank);
  
  SYS_T::print_fatal_if(size != PartBasic->get_cpu_size(), "Error: number of processors does not match with prepost! \n");

  SYS_T::commPrint("===> %d processor(s) are assigned for:", size);

  // Allocate the quadrature rule and element container
  IQuadPts * quad = nullptr;
  FEAElement * element = nullptr; 

  if( GMIptr->get_elemType() == 501 )
  {
    quad = new QuadPts_vis_tet4();
    element = new FEAElement_Tet4( quad-> get_num_quadPts() );
  }
  else if( GMIptr->get_elemType() == 502 )
  {
    quad = new QuadPts_vis_tet10_v2();
    element = new FEAElement_Tet10_v2( quad-> get_num_quadPts() );
  }
  else SYS_T::print_fatal( "Error: unsupported element type \n" );

  // Print the sampling points on screen
  quad -> print_info();

  // Create the visualization data object
  IVisDataPrep * visprep = new VisDataPrep_NS();

  visprep->print_info();
 
  // Allocate the container to store the solution values into physical fields.
  double ** solArrays = new double * [visprep->get_ptarray_size()];
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    solArrays[ii] = new double [pNode->get_nlocghonode() * visprep->get_ptarray_comp_length(ii)];

  // VTK writer 
  VTK_Writer_NS * vtk_w = new VTK_Writer_NS( GMIptr->get_nElem(), 
      GMIptr->get_nLocBas(), element_part_file );

  std::ostringstream time_index;

  // ==== WOMERSLEY CHANGES BEGIN ====
  const int nElem   = locElem->get_nlocalele();
  const int nLocBas = GMIptr->get_nLocBas();

  // Initialize velo and pres errors in L2, H1 norms
  double err_velo_l2 = 0.0, err_velo_h1 = 0.0, err_pres_l2 = 0.0, err_pres_h1 = 0.0;

  // Allocate element arrays 
  int * IEN_e = new int [nLocBas];
  double * R  = new double [nLocBas];
  double * Rx = new double [nLocBas];
  double * Ry = new double [nLocBas];
  double * Rz = new double [nLocBas];

  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];

  double * loc_p = new double [nLocBas];
  double * loc_u = new double [nLocBas];
  double * loc_v = new double [nLocBas];
  double * loc_w = new double [nLocBas];
  // ==== WOMERSLEY CHANGES END ====

  for(int time = time_start; time<=time_end; time+= time_step)
  {
    std::string name_to_read(sol_bname);
    std::string name_to_write(out_bname);
    time_index.str("");
    time_index<< 900000000 + time;
    name_to_read.append(time_index.str());
    name_to_write.append(time_index.str());

    SYS_T::commPrint("Time %d: Read %s and Write %s \n",
        time, name_to_read.c_str(), name_to_write.c_str() );

    visprep->get_pointArray(name_to_read, anode_mapping_file, pnode_mapping_file,
        pNode, GMIptr, dof, solArrays);

    
    vtk_w->writeOutput( fNode, locIEN, locElem,
        visprep, element, quad, solArrays,
        rank, size, time * dt, sol_bname, out_bname, name_to_write, isXML );

    // ==== WOMERSLEY CHANGES BEGIN ====
    if(time == manu_sol_time)
    {
      PostVectSolution * pvsol = new PostVectSolution(name_to_read,
          anode_mapping_file, pnode_mapping_file, pNode, GMIptr, dof);

      for(int ee=0; ee<nElem; ++ee)
      {
        // Build element
        locIEN  -> get_LIEN_e(ee, IEN_e);
        fNode   -> get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);
        element -> buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );

        // Obtain element solution vector
        pvsol -> get_esol(0, nLocBas, IEN_e, loc_p);
        pvsol -> get_esol(1, nLocBas, IEN_e, loc_u);
        pvsol -> get_esol(2, nLocBas, IEN_e, loc_v);
        pvsol -> get_esol(3, nLocBas, IEN_e, loc_w);

        // Compute errors
        err_velo_l2 += POST_T_NS::get_velo_l2_error( loc_u, loc_v, loc_w,
            element, ectrl_x, ectrl_y, ectrl_z, quad, R, time * dt);

        err_velo_h1 += POST_T_NS::get_velo_h1_error( loc_u, loc_v, loc_w,
            element, ectrl_x, ectrl_y, ectrl_z, quad, R, Rx, Ry, Rz, time * dt);

        err_pres_l2 += POST_T_NS::get_pres_l2_error( loc_p,
            element, ectrl_x, ectrl_y, ectrl_z, quad, R, time * dt );

        err_pres_h1 += POST_T_NS::get_pres_h1_error( loc_p,
            element, ectrl_x, ectrl_y, ectrl_z, quad, R, Rx, Ry, Rz, time * dt );
      }

      delete pvsol;
    }
    // ==== WOMERSLEY CHANGES END ====
  }

  // ==== WOMERSLEY CHANGES BEGIN ====
  MPI_Barrier(PETSC_COMM_WORLD);
  double velo_l2 = 0.0, velo_h1 = 0.0, pres_l2 = 0.0, pres_h1 = 0.0;
  MPI_Reduce(&err_velo_l2, &velo_l2, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&err_velo_h1, &velo_h1, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&err_pres_l2, &pres_l2, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
  MPI_Reduce(&err_pres_h1, &pres_h1, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);

  velo_l2 = sqrt(velo_l2); velo_h1 = sqrt(velo_h1 + velo_l2*velo_l2);
  pres_l2 = sqrt(pres_l2); pres_h1 = sqrt(pres_h1 + pres_l2*pres_l2);

  PetscPrintf(PETSC_COMM_WORLD, "Abs Error in L2 norm of velo at t=%f is : %e \n", manu_sol_time * dt, velo_l2);
  PetscPrintf(PETSC_COMM_WORLD, "Abs Error in H1 norm of velo at t=%f is : %e \n", manu_sol_time * dt, velo_h1);

  PetscPrintf(PETSC_COMM_WORLD, "Abs Error in L2 norm of pres at t=%f is : %e \n", manu_sol_time * dt, pres_l2);
  PetscPrintf(PETSC_COMM_WORLD, "Abs Error in H1 norm of pres at t=%f is : %e \n", manu_sol_time * dt, pres_h1);
  // ==== WOMERSLEY CHANGES END ====

  // ===== Clean the memory =====
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    delete [] solArrays[ii];
  delete [] solArrays;
  
  delete fNode; delete locIEN; delete GMIptr; delete PartBasic; delete locElem;
  delete pNode; delete quad; delete element; delete visprep; delete vtk_w;

  // ==== WOMERSLEY CHANGES BEGIN ====
  delete [] IEN_e; delete [] R; delete [] Rx; delete [] Ry; delete [] Rz;
  delete [] ectrl_x; delete [] ectrl_y; delete [] ectrl_z;
  delete [] loc_p; delete [] loc_u; delete [] loc_v; delete [] loc_w;
  // ==== WOMERSLEY CHANGES END ====

  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
