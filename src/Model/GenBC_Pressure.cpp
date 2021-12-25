#include "GenBC_Pressure.hpp"

GenBC_Pressure::GenBC_Pressure( const std::string &lpn_filename, const double &init_time )
{
  num_ebc = 0;
  coef_a.clear(); coef_b.clear(); ebc_ids.clear();
  num_of_mode.clear(); w.clear(); period.clear();

  if( SYS_T::file_exist( lpn_filename ) )
  {
    std::ifstream reader;
    reader.open( lpn_filename.c_str(), std::ifstream::in );

    std::istringstream sstrm;
    std::string sline;
    std::string bc_type;

    // The first non-commented, non-empty line should be
    // Pressure num_ebc
    while( std::getline(reader, sline) )
    {
      if( sline[0] !='#' && !sline.empty() )
      {
        sstrm.str(sline);
        sstrm >> bc_type >> num_ebc;
        sstrm.clear();
        break;
      }
    }

    if( bc_type.compare("Pressure") != 0 && bc_type.compare("PRESSURE") != 0 )
      SYS_T::print_fatal( "GenBC_Pressure Error: BC type in %s should be Pressure.\n", lpn_filename.c_str() );

    coef_a.resize(num_ebc); coef_b.resize(num_ebc);
    num_of_mode.resize(num_ebc); w.resize(num_ebc); period.resize(num_ebc);

    // Read in num_of_mode, w, period, coef_a, and coef_b per ebc
    for(int ii=0; ii<num_ebc; ++ii)
    {
      while( std::getline(reader, sline) )
      {
        // face_id num_of_mode w period
        if( sline[0] !='#' && !sline.empty() )
        {
          sstrm.str(sline);
          int face_id;

          sstrm >> face_id >> num_of_mode[ii] >> w[ii] >> period[ii];
          sstrm.clear();

          ebc_ids.push_back( face_id );
          break;
        }
      }

      // Check the compatibility of period and w. If the difference
      // is larger than 0.01, print a warning message
      if( std::abs(2.0 * MATH_T::PI / period[ii] - w[ii] ) >= 0.01 ) SYS_T::commPrint( "\nGenBC_Pressure WARNING: ebc_id %d incompatible period and w, \n2xpi/period = %e and w = %e.\n", ebc_ids[ii], 2.0*MATH_T::PI/period[ii], w[ii] );

      coef_a[ii].clear(); coef_b[ii].clear();

      while( std::getline(reader, sline) )
      {
        // coef_a
        if( sline[0] !='#' && !sline.empty() )
        {
          sstrm.str(sline);
          double temp_coef;
          while( sstrm >> temp_coef ) coef_a[ii].push_back( temp_coef );

          sstrm.clear();
          break;
        }
      }

      VEC_T::shrink2fit( coef_a[ii] );

      if( static_cast<int>(coef_a[ii].size()) != num_of_mode[ii]+1 ) SYS_T::print_fatal( "GenBC_Pressure Error: ebc_id %d a-coefficients in %s incompatible with the given number of modes.\n", ebc_ids[ii], lpn_filename.c_str() );

      while( std::getline(reader, sline) )
      {
        // coef_b
        if( sline[0] !='#' && !sline.empty() )
        {
          sstrm.str(sline);
          double temp_coef;
          while( sstrm >> temp_coef ) coef_b[ii].push_back( temp_coef );

          sstrm.clear();
          break;
        }
      }

      VEC_T::shrink2fit( coef_b[ii] );

      if( static_cast<int>(coef_b[ii].size()) != num_of_mode[ii]+1 ) SYS_T::print_fatal( "GenBC_Pressure Error: ebc_id %d b-coefficients in %s incompatible with the given number of modes.\n", ebc_ids[ii], lpn_filename.c_str() );
    }

    if( (int) ebc_ids.size() != num_ebc ) SYS_T::print_fatal("Error: GenBC_Pressure the input file %s does not contain complete data for outlet faces. \n", lpn_filename.c_str());

    // Finish reading the file and close it
    reader.close();

    // Initialize P0 by setting it to be get_P at time = 0.0
    P0.resize( num_ebc );
    for(int ii = 0; ii < num_ebc; ++ii)
      P0[ii] = get_P( ii, 0.0, 0.0, init_time );
  }
}

GenBC_Pressure::~GenBC_Pressure()
{
  VEC_T::clean(coef_a); VEC_T::clean(coef_b);
  VEC_T::clean(num_of_mode); VEC_T::clean(w); VEC_T::clean(period);
}

void GenBC_Pressure::print_info() const
{
  SYS_T::commPrint("===> GenBC_Pressure: \n");

  for(int ii=0; ii<num_ebc; ++ii)
  {
    SYS_T::commPrint("     ebc_id = %d", ebc_ids[ii]);
    SYS_T::commPrint(" w = %e, period =%e \n", w[ii], period[ii]);
    SYS_T::commPrint("     a[0] + Sum{ a[i] cos(i x w x t) + b[i] sin(i x w x t) }, for i = 1,...,%d. \n", num_of_mode[ii]);
    for(int jj=0; jj<=num_of_mode[ii]; ++jj)
      SYS_T::commPrint("     i = %d, a = %e, b = %e \n", jj, coef_a[ii][jj], coef_b[ii][jj]);
  }
}

double GenBC_Pressure::get_P( const int &idx, const double &dot_Q, const double &Q,
       const double &time ) const
{
  double PP = coef_a[idx][0];
  for( int ii = 1; ii <= num_of_mode[idx]; ++ii )
    PP += coef_a[idx][ii] * cos( ii*w[idx]*time ) + coef_b[idx][ii] * sin( ii*w[idx]*time );

  return PP;
}

// EOF
