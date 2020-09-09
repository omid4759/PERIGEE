// ==================================================================
// post_compare_manu.cpp
// ------------------------------------------------------------------
// This is a postprocessing driver for computing the solution error
// given a manufactured solution.
//
// Date: Sep. 8 2020
// ==================================================================

#include "Tet_Tools.hpp"
#include "Post_error_ns.hpp"

int main( int argc, char * argv[] )
{
  std::string geo_file, wall_file;    // volumetric and wall files
  int elemType;                       // 501 for linear tet; 502 for quadratic tet
  int dof;
  double fluid_mu;

  // Enforce serial execution
  PetscMPIInt size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  SYS_T::print_fatal_if(size!=1, "ERROR: post_compare_manu is a serial program! \n");

  // Read preprocessor_cmd.h5
  hid_t prepcmd_file = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );
  cmd_h5r -> read_string("/", "geo_file",      geo_file);
  cmd_h5r -> read_string("/", "sur_file_wall", wall_file);
  elemType = cmd_h5r -> read_intScalar   ("/", "elemType");
  dof      = cmd_h5r -> read_doubleScalar("/", "dofNum");
  delete cmd_h5r; H5Fclose(prepcmd_file);

  // Read solver_cmd.h5
  prepcmd_file = H5Fopen("solver_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  cmd_h5r = new HDF5_Reader( prepcmd_file );
  fluid_mu = cmd_h5r -> read_doubleScalar("/", "fl_mu");
  delete cmd_h5r; H5Fclose(prepcmd_file);

  std::string sol_name;                // solution filename
  double sol_time;                     // time for error calculation

  // ===== Command Line Arguments ===== 
  SYS_T::GetOptionReal(  "-sol_time", sol_time);
  SYS_T::GetOptionString("-sol_name", sol_name);

  cout << "==== Command Line Arguments ====" << endl;
  cout << " sol_time: "  << sol_time  << endl;
  cout << " sol_name: "  << sol_name  << endl;
  cout << "----------------------------------\n";
  cout << " geo_file: "  << geo_file  << endl;
  cout << " wall_file: " << wall_file << endl;
  cout << " elemType: "  << elemType  << endl;
  cout << " fl_mu: "     << fluid_mu  << endl;
  cout <<"==== Command Line Arguments ===="<<endl;

  // Make sure the files exist on disk
  SYS_T::file_check( geo_file.c_str() );
  SYS_T::file_check( wall_file.c_str() );

  // Make sure the files exist on disk
  SYS_T::file_check( geo_file.c_str() );
  SYS_T::file_check( wall_file.c_str() );

  // Read in the volumetric mesh info
  int v_nFunc, v_nElem;
  std::vector<int> v_vecIEN;
  std::vector<double> v_ctrlPts;

  TET_T::read_vtu_grid(geo_file.c_str(), v_nFunc, v_nElem, v_ctrlPts, v_vecIEN);
  cout<<endl<<"Volumetric mesh contains "<<v_nElem<<" elements and "<<v_nFunc<<" vertices.\n";

  // Read in the wall surface mesh info
  int nFunc, nElem;
  std::vector<double> ctrlPts;
  std::vector<int> vecIEN, global_node_idx, global_ele_idx;

  if(elemType == 501)
  {
    TET_T::read_vtp_grid( wall_file.c_str(), nFunc, nElem, ctrlPts, vecIEN,
      global_node_idx, global_ele_idx );    
  }
  else if(elemType == 502)
  {
    TET_T::read_vtu_grid( wall_file.c_str(), nFunc, nElem, ctrlPts, vecIEN,
      global_node_idx, global_ele_idx );
  }
  else SYS_T::print_fatal("Error: unknown elem type.\n");
  cout<<"Wall mesh contains "<<nElem<<" elements and "<<nFunc<<" vertices.\n";



}

// EOF
