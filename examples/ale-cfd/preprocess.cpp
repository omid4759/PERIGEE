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

using std::cout;
using std::endl;

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
  const std::string yaml_file("ns_preprocess.yml");

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

  return EXIT_SUCCESS;
}

// EOF
