// ============================================================================
// This is a preprocessor code for ALE-CFD analysis
// 
// Date Created: Dec. 01 2023
// ============================================================================
#include "Math_Tools.hpp"
#include "Mesh_Tet.hpp"
#include "Mesh_FEM.hpp"
#include "IEN_FEM.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_FEM.hpp"
#include "NodalBC.hpp"

#include "yaml-cpp/yaml.h"

int main( int argc, char * argv[] )
{
  // Set number of threads and  print info of OpenMP
  SYS_T::print_omp_info();
  SYS_T::set_omp_num_threads();

  // Clean the potentially pre-existing hdf5 files in the job folder
  SYS_T::execute("rm -rf part_p*.h5");
  SYS_T::execute("rm -rf preprocessor_cmd.h5");

  // Define basic problem settins
  constexpr int dofNum = 7; // degree-of-freedom for the physical problem
  constexpr int dofMat = 4; // degree-of-freedom in the matrix problem

  // Yaml options
  const std::string yaml_file("ale_cfd_preprocess.yml");

  // Check if the yaml file exist on disk
  SYS_T::file_check(yaml_file);

  YAML::Node paras = YAML::LoadFile( yaml_file );

  const int elemType                  = paras["elem_type"].as<int>();
  const int num_inlet                 = paras["num_inlet"].as<int>();
  const int num_outlet                = paras["num_outlet"].as<int>();
  const std::string geo_file          = paras["geo_file"].as<std::string>();
  const std::string sur_file_in_base  = paras["sur_file_in_base"].as<std::string>();
  const std::string sur_file_wall     = paras["sur_file_wall"].as<std::string>();
  const std::string sur_file_out_base = paras["sur_file_out_base"].as<std::string>();
  const std::string part_file         = paras["part_file"].as<std::string>();
  const int cpu_size                  = paras["cpu_size"].as<int>();
  const int in_ncommon                = paras["in_ncommon"].as<int>();
  const bool isDualGraph              = paras["is_dualgraph"].as<bool>();

  if( elemType != 501 && elemType != 502 && elemType != 601 && elemType != 602 ) SYS_T::print_fatal("ERROR: unknown element type %d.\n", elemType);

  // Print the command line arguments
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -elem_type: "<<elemType<<endl;
  cout<<" -num_outlet: "<<num_outlet<<endl;
  cout<<" -geo_file: "<<geo_file<<endl;
  cout<<" -sur_file_in_base: "<<sur_file_in_base<<endl;
  cout<<" -sur_file_wall: "<<sur_file_wall<<endl;
  cout<<" -sur_file_out_base: "<<sur_file_out_base<<endl;
  cout<<" -part_file: "<<part_file<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph) cout<<" -isDualGraph: true \n";
  else cout<<" -isDualGraph: false \n";
  cout<<"---- Problem definition ----\n";
  cout<<" dofNum: "<<dofNum<<endl;
  cout<<" dofMat: "<<dofMat<<endl;
  cout<<"====  Command Line Arguments/ ===="<<endl;

  // Check if the geometry files exist on disk
  SYS_T::file_check(geo_file); cout<<geo_file<<" found. \n";

  SYS_T::file_check(sur_file_wall); cout<<sur_file_wall<<" found. \n";

  // Generate the inlet file names and check existance
  std::vector< std::string > sur_file_in;
  sur_file_in.resize( num_inlet );

  for(int ii=0; ii<num_inlet; ++ii)
  {
    if(elemType == 501 || elemType == 601)
      sur_file_in[ii] = SYS_T::gen_capfile_name( sur_file_in_base, ii, ".vtp" );
    else if(elemType == 502 || elemType == 602)
      sur_file_in[ii] = SYS_T::gen_capfile_name( sur_file_in_base, ii, ".vtu" );
    else
      SYS_T::print_fatal("Error: unknown element type occurs when generating the inlet file names. \n");

    SYS_T::file_check(sur_file_in[ii]);
    cout<<sur_file_in[ii]<<" found. \n";
  }

  // Generate the outlet file names and check existance
  std::vector< std::string > sur_file_out;
  sur_file_out.resize( num_outlet );

  for(int ii=0; ii<num_outlet; ++ii)
  {
    if(elemType == 501 || elemType == 601)
      sur_file_out[ii] = SYS_T::gen_capfile_name( sur_file_out_base, ii, ".vtp" );
    else if(elemType == 502 || elemType == 602)
      sur_file_out[ii] = SYS_T::gen_capfile_name( sur_file_out_base, ii, ".vtu" );
    else
      SYS_T::print_fatal("Error: unknown element type occurs when generating the outlet file names. \n");

    SYS_T::file_check(sur_file_out[ii]);
    cout<<sur_file_out[ii]<<" found. \n";
  }

  // Record the problem setting into a HDF5 file: preprocessor_cmd.h5
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("num_inlet", num_inlet);
  cmdh5w->write_intScalar("num_outlet", num_outlet);
  cmdh5w->write_intScalar("cpu_size", cpu_size);
  cmdh5w->write_intScalar("in_ncommon", in_ncommon);
  cmdh5w->write_intScalar("dofNum", dofNum);
  cmdh5w->write_intScalar("dofMat", dofMat);
  cmdh5w->write_intScalar("elemType", elemType);
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("sur_file_in_base", sur_file_in_base);
  cmdh5w->write_string("sur_file_out_base", sur_file_out_base);
  cmdh5w->write_string("sur_file_wall", sur_file_wall);
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);





  return EXIT_SUCCESS;
}

// EOF
