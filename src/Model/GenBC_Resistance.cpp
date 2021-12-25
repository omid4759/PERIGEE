#include "GenBC_Resistance.hpp"

GenBC_Resistance::GenBC_Resistance( const std::string &lpn_filename )
{
  num_ebc = 0;
  time.clear(); resis.clear(); pres_offset.clear();
  Q0.clear(); P0.clear(); ebc_ids.clear();

  // Now read the values of resis_R and resis_Pd from disk file lpn_filename
  if( SYS_T::file_exist( lpn_filename ) )
  {
    std::ifstream reader;
    reader.open( lpn_filename.c_str(), std::ifstream::in );

    std::istringstream sstrm;
    std::string sline;

    std::string bc_type;

    // The first non-commented line should contain
    // bc_type [Resistance, RCR] number-of-ebc
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
    if( bc_type.compare("Resistance") != 0 
        && bc_type.compare("resistance") != 0 
        && bc_type.compare("RESISTANCE") != 0 )
    {
      SYS_T::print_fatal("Error: the outflow model in %s does not match GenBC_Resistance.\n", lpn_filename.c_str());
    }

    time.resize(num_ebc); resis.resize(num_ebc); pres_offset.resize(num_ebc);
    Q0.resize(num_ebc); P0.resize(num_ebc);
 
    // Read files for each ebc to define the parameters for LPN
    for( int ii=0; ii<num_ebc; ++ii )
    {
      while( std::getline(reader, sline) )
      {
        if( sline[0] != '#' && !sline.empty() )
        {
          sstrm.str( sline );
          int face_id, num_timept;
          sstrm >> face_id >> num_timept;

          SYS_T::print_fatal_if( num_timept < 2,
              "Error: GenBC_Resistance number of timepoints must be >= 2. \n" );

          time[ii].resize( num_timept );
          resis[ii].resize( num_timept );
          pres_offset[ii].resize( num_timept );

          for( int tt=0; tt<num_timept; ++tt)
          {
            std::getline(reader, sline);
            sstrm.str( sline );
            sstrm >> time[ii][tt] >> resis[ii][tt] >> pres_offset[ii][tt];
            sstrm.clear();
          }

          sstrm.clear();
          ebc_ids.push_back( face_id );
          break;
        }
      }
    }

    if( (int) ebc_ids.size() != num_ebc ) SYS_T::print_fatal("Error: GenBC_Resistance the input file %s does not contain complete data for outlet faces. \n", lpn_filename.c_str());

    reader.close();

    SYS_T::commPrint( "===> GenBC_Resistance data are read in from %s.\n", lpn_filename.c_str() );

    // Set zero initial value.
    for(int ii=0; ii<num_ebc; ++ii)
    {
      Q0[ii] = 0.0;
      P0[ii] = 0.0;
    }
  }
}

GenBC_Resistance::~GenBC_Resistance()
{}

void GenBC_Resistance::print_info() const
{
  SYS_T::commPrint("===> GenBC_Resistance : \n");
  for(int ii=0; ii<num_ebc; ++ii)
  {
    SYS_T::commPrint( "     ebcid = %d, num_timept = %d \n", ebc_ids[ii], (int) resis[ii].size() );
    if( resis[ii].size() < 10)
    {
      for( int tt = 0; tt < (int) resis[ii].size(); ++tt )
        SYS_T::commPrint( "     t = %e, R = %e, Pd = %e \n", time[ii][tt], resis[ii][tt], pres_offset[ii][tt] );
    }
  }
}

double GenBC_Resistance::get_m( const int &ii, const double &dot_Q, const double &Q,
    const double &curr_time ) const
{
  const int num_timept = time[ii].size();
  const double period = time[ii][num_timept - 1];
  const double time_fmod = std::fmod( curr_time, period );

  const int tt_n = get_interp_ttn( ii, curr_time );
  const double dt = time[ii][tt_n + 1] - time[ii][tt_n]; 

  const double res = ( resis[ii][tt_n] * ( time[ii][tt_n + 1] - time_fmod ) +
      resis[ii][tt_n + 1] * ( time_fmod - time[ii][tt_n] ) ) / dt;

  return res;
}

double GenBC_Resistance::get_P( const int &ii, const double &dot_Q, const double &Q,
       const double &curr_time ) const
{
  const int num_timept = time[ii].size();
  const double period = time[ii][num_timept - 1];
  const double time_fmod = std::fmod( curr_time, period );

  const int tt_n = get_interp_ttn( ii, curr_time );
  const double dt = time[ii][tt_n + 1] - time[ii][tt_n]; 

  const double res = ( resis[ii][tt_n] * ( time[ii][tt_n + 1] - time_fmod ) +
      resis[ii][tt_n + 1] * ( time_fmod - time[ii][tt_n] ) ) / dt;

  const double pd = ( pres_offset[ii][tt_n] * ( time[ii][tt_n + 1] - time_fmod ) +
      pres_offset[ii][tt_n + 1] * ( time_fmod - time[ii][tt_n] ) ) / dt;

  return res * Q + pd;
}


int GenBC_Resistance::get_interp_ttn( const int &ii, const double &curr_time ) const
{
  const int num_timept = time[ii].size();
  const double period = time[ii][num_timept - 1];
  const double time_fmod = std::fmod( curr_time, period );

  for(int tt = 0; tt < num_timept - 1; ++tt)
  {
    if( time_fmod >= time[ii][tt] && time_fmod < time[ii][tt + 1] )
      return tt;
  }

  SYS_T::print_fatal("Error: GenBC_Resistance interpolation tt_n not found.\n");
  return -1;
}

// EOF
