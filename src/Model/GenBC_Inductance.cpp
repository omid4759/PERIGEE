#include "GenBC_Inductance.hpp"

GenBC_Inductance::GenBC_Inductance( const std::string &lpn_filename )
{
  num_ebc = 0;
  induct.clear(); pres_offset.clear(); ebc_ids.clear();
  Q0.clear(); P0.clear();

  // Now read the values of induct_L and induct_p from disk file lpn_filename
  if( SYS_T::file_exist( lpn_filename ) )
  {
    std::ifstream reader;
    reader.open( lpn_filename.c_str(), std::ifstream::in );

    std::istringstream sstrm;
    std::string sline;

    std::string bc_type;

    // The first non-commented line should contain
    // bc_type [Resistance, RCR, Inductance] number-of-ebc
    while( std::getline(reader, sline) )
    {
      if( sline[0] != '#' && !sline.empty() )
      {
        sstrm.str(sline);
        sstrm >> bc_type >> num_ebc;
        sstrm.clear();
        break;
      }
    }

    // Allocate the distal pressure and resistance vector
    if( bc_type.compare("Inductance") != 0 
        && bc_type.compare("inductance") != 0 
        && bc_type.compare("INDUCTANCE") != 0 )
    {
      SYS_T::print_fatal("Error: the outflow model in %s does not match GenBC_Inductance.\n", lpn_filename.c_str());
    }
 
    Q0.resize( num_ebc ); P0.resize( num_ebc );

    // Read files for each ebc to define the parameters for LPN
    while( std::getline(reader, sline) )
    {
      if( sline[0] != '#' && !sline.empty() )
      {
        sstrm.str( sline );
        int face_id;
        double L, pd;

        sstrm >> face_id >> L >> pd;
        sstrm.clear();

        induct.push_back( L );
        pres_offset.push_back( pd );
        ebc_ids.push_back( face_id );  
      }
    }

    if( (int) ebc_ids.size() != num_ebc ) SYS_T::print_fatal("Error: GenBC_Inductance the input file %s does not contain complete data for outlet faces. \n", lpn_filename.c_str());

    reader.close();

    SYS_T::commPrint( "===> GenBC_Inductance data are read in from %s.\n", lpn_filename.c_str() );

    // Set zero initial value.
    for(int ii=0; ii<num_ebc; ++ii)
    {
      Q0[ii] = 0.0;
      P0[ii] = 0.0;
    }
  }
}


GenBC_Inductance::~GenBC_Inductance()
{}


void GenBC_Inductance::print_info() const
{
  SYS_T::commPrint("===> GenBC_Inductance : \n");
  for(int ii=0; ii<num_ebc; ++ii)
    SYS_T::commPrint( "     ebcid = %d, L = %e, p = %e \n", ebc_ids[ii], induct[ii], pres_offset[ii] );
}

// EOF
