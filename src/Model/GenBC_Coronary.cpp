#include "GenBC_Coronary.hpp"

GenBC_Coronary::GenBC_Coronary( const std::string &lpn_filename,
    const int &in_N, const double &dt3d, const int &in_index,
    const std::string &in_lpn_sol_file )
: num_odes(2), N( in_N ), h( dt3d/static_cast<double>(N) ),
  lpn_sol_file( in_lpn_sol_file ), absTol( 1.0e-8 ), relTol( 1.0e-5 )
{
  num_ebc = 0;
  Ra.clear(); Ca.clear(); Ra_micro.clear();
  Cim.clear(); Rv.clear(); Pd.clear();
  num_Pim_data.clear(); alpha_Pim.clear(); ebc_ids.clear();
  Q0.clear(); Pi0.clear(); Time_data.clear(); Pim_data.clear();
  der_Pim_data.clear(); prev_0D_sol.clear(); restart_0D_sol.clear();
  dPimdt_k1.clear(); dPimdt_k2.clear(); dPimdt_k3.clear();

  // Now read the lpn input file for num_ebc and coronary model
  // parameters (Ra, Ca, Ra_micro, Cim, Rv, Pd and Pim)
  if( SYS_T::file_exist( lpn_filename ) )
  {
    std::ifstream reader;
    reader.open( lpn_filename.c_str(), std::ifstream::in );

    std::istringstream sstrm;
    std::string sline;
    std::string bc_type;

    // The first non-commented line should be Coronary num_ebc
    while( std::getline(reader, sline) )
    {
      if( sline[0] != '#' && !sline.empty() )
      {
        sstrm.str(sline);
        sstrm >> bc_type;
        sstrm >> num_ebc;
        sstrm.clear();
        break;
      }
    }

    // Check the file's bc_type matches Coronary
    if( bc_type.compare("Coronary") != 0
        && bc_type.compare("CORONARY") != 0
        && bc_type.compare("coronary") != 0 )
    {
      SYS_T::print_fatal("Error: the outflow model in %s does not match GenBC_Coronary.\n", lpn_filename.c_str());
    }

    // Allocate the containers for model parameters
    Q0.resize( num_ebc );             Pi0.resize( num_ebc );
    Time_data.resize( num_ebc );      Pim_data.resize( num_ebc );
    der_Pim_data.resize( num_ebc );   prev_0D_sol.resize( num_ebc );
    restart_0D_sol.resize( num_ebc ); dPimdt_k1.resize( num_ebc );
    dPimdt_k2.resize( num_ebc );      dPimdt_k3.resize( num_ebc );

    for( int ii =0; ii<num_ebc; ++ii )
    {
      prev_0D_sol[ii].resize( num_odes );
      restart_0D_sol[ii].resize( num_odes );
      Pi0[ii].resize( num_odes );

      dPimdt_k1[ii].resize( N+1 );
      dPimdt_k2[ii].resize( N );
      dPimdt_k3[ii].resize( N );
    }

    // Read files for each ebc to set the values of Ra, Ca, Ra_micro,
    // Cim, Rv, Pd, num_Pim_data, and alpha_Pim
    for( int ii=0; ii<num_ebc; ++ii ) 
    {
      while( std::getline(reader, sline) )
      {
        if( sline[0] != '#' && !sline.empty() )
        {
          sstrm.str( sline );
          int face_id, numpimdata;
          double ra, ca, ramicro, cim, rv, pd, alphapim;

          sstrm >> face_id >> ra >> ca >> ramicro >> cim >> rv >> pd >> numpimdata >> alphapim;

          Ra.push_back( ra );
          Ca.push_back( ca );
          Ra_micro.push_back( ramicro );
          Cim.push_back( cim );
          Rv.push_back( rv );
          Pd.push_back( pd );
          num_Pim_data.push_back( numpimdata );
          alpha_Pim.push_back( alphapim );

          SYS_T::print_fatal_if( numpimdata < 2,
              "Error: number of Pim data must be >= 2 for coronary BC. \n" );

          // Resize the IntraMyocardial data
          Time_data[ii].resize( numpimdata );
          Pim_data[ii].resize(  numpimdata );
          der_Pim_data[ii].resize( numpimdata );

          sstrm.clear();

          for(int tt = 0; tt < numpimdata; ++tt)
          {
            std::getline(reader, sline);
            sstrm.str( sline );
            sstrm >> Time_data[ii][tt];
            sstrm >> Pim_data[ii][tt];

            Pim_data[ii][tt] = Pim_data[ii][tt] * alpha_Pim[ii];
            sstrm.clear();
          }

          GENBC_T::set_pchip( num_Pim_data[ii], Time_data[ii], Pim_data[ii], der_Pim_data[ii] );
          get_dPim_dt( ii, 0.0, N * h );

          SYS_T::print_fatal_if( Time_data[ii][0] > 0.0, "Error: Pim data does not start from 0.\n" );

          sstrm.clear();
          ebc_ids.push_back( face_id );
          break;
        }
      }
    }

    SYS_T::print_fatal_if( (int) ebc_ids.size() != num_ebc, "Error: GenBC_Coronary the input file %s does not contain complete data for outlet faces.\n", lpn_filename.c_str() );

    reader.close();

    SYS_T::commPrint( "===> GenBC_Coronary data are read in from %s.\n", lpn_filename.c_str() );

    // Set zero initial values. These will be reset based on the 3D solutions.
    for(int ii=0; ii<num_ebc; ++ii)
    {
      Q0[ii] = 0.0;
      for(int jj=0; jj<num_odes; ++jj)
      {
        Pi0[ii][jj] = 0.0;
        prev_0D_sol[ii][jj] = 0.0;
        restart_0D_sol[ii][jj] = 0.0;
      }

      // Make sure C and R are nonzero
      SYS_T::print_fatal_if( Ca[ii] == 0.0, "Error: GenBC_Coronary Ca cannot be zero.\n" );
      SYS_T::print_fatal_if( Cim[ii] == 0.0, "Error: GenBC_Coronary Cim cannot be zero.\n" );
      SYS_T::print_fatal_if( Ra_micro[ii] == 0.0, "Error: GenBC_Coronary Ra_micro cannot be zero.\n" );
      SYS_T::print_fatal_if( Rv[ii] == 0.0, "Error: GenBC_Coronary Rv cannot be zero.\n" );
    }

    if( in_index == 0 )
    {
      std::ofstream ofile;
      ofile.open( lpn_sol_file.c_str(), std::ofstream::out | std::ofstream::trunc );
      ofile<<"# Time index"<<'\t'<<"Time"<<'\t'<<"0D_sol[i][j] 0D_sol[i][j+1] ..."<<'\t'<<"i=0...num_ebc-1, j=0...num_odes-1"<<'\n';
      ofile.close();
    }
    // if restart index > 0, read initial 0D solutions from coronary_sol.txt
    else
    {
      reader.open( lpn_sol_file.c_str(), std::ifstream::in );
      bool is_0D_sol_found = false;

      while( std::getline(reader, sline) )
      {
        if( sline[0] != '#' && !sline.empty() )
        {
          sstrm.str(sline);
          int temp_index;
          double temp_time;
          sstrm >> temp_index;
          sstrm >> temp_time;
          if( temp_index == in_index )
          {
            is_0D_sol_found = true;
            for( int ii=0; ii<num_ebc; ++ii )
            {
              for( int jj=0; jj<num_odes; ++jj )
              {
                sstrm >> restart_0D_sol[ii][jj];
                Pi0[ii][jj] = restart_0D_sol[ii][jj];
 
                // Precalculate dPimdt values needed for integrating coronary ODEs.
                get_dPim_dt(ii, temp_time, temp_time + N * h);
              }
            }
            break;
          }
        }
      }

      reader.close();
      SYS_T::print_fatal_if(!is_0D_sol_found, "Error: GenBC_Coronary cannot find 0D solutions for restart index = %d .\n", in_index);
    }
  }
}

GenBC_Coronary::~GenBC_Coronary()
{}

void GenBC_Coronary::print_info() const
{
  SYS_T::commPrint( "===> Coronary model: N = %d, h = %e, num_ebc = %d \n", N, h, num_ebc );

  for(int ii=0; ii<num_ebc; ++ii)
    SYS_T::commPrint( "     ebcid = %d, Ra = %e, Ca = %e, Ra_micro = %e,Rim = %e, Rv = %e, Pd = %e \n",
        ebc_ids[ii], Ra[ii], Ca[ii],Ra_micro[ii],Cim[ii], Rv[ii], Pd[ii] );
}

void GenBC_Coronary::write_0D_sol( const int &curr_index, const double &curr_time ) const
{
  std::ofstream ofile;
  ofile.open( lpn_sol_file.c_str(), std::ofstream::out | std::ofstream::app );

  ofile << curr_index << '\t' << curr_time;

  for( int ii=0; ii<num_ebc; ++ii )
    for( int jj=0; jj<num_odes; ++jj ) ofile << '\t' << std::scientific << std::setprecision(16) << std::setw(20) << prev_0D_sol[ii][jj];

  ofile << '\n';
  ofile.close();
}

double GenBC_Coronary::get_m( const int &ii, const double &in_dot_Q, const double &in_Q ) const
{
  double diff = std::abs(in_Q) * relTol;

  if( diff < absTol ) diff = absTol;

  const double left  = get_P(ii, in_dot_Q, in_Q + 0.5 * diff);
  const double right = get_P(ii, in_dot_Q, in_Q - 0.5 * diff);

  return (left - right) / diff;
}

double GenBC_Coronary::get_P( const int &ii, const double &in_dot_Q,
    const double &in_Q, const double &time ) const
{
  double output_P = 0.0;

  const double fac13 = 1.0 / 3.0;
  const double fac23 = 2.0 / 3.0;
  const double fac18 = 1.0 / 8.0;
  const double fac38 = 3.0 / 8.0;

  // Each coronary face is governed by num_odes = 2 ODEs.
  // initial pressures at Ca and Cim
  double pi_m[num_odes];

  // auxiliary variables for RK4
  double K1[num_odes], K2[num_odes], K3[num_odes], K4[num_odes];
  double pi_tmp[num_odes];

  for(int jj=0; jj<num_odes; ++jj)
  {
    pi_m[jj] = Pi0[ii][jj];
    K1[jj] = 0.0;
    K2[jj] = 0.0;
    K3[jj] = 0.0;
    K4[jj] = 0.0;
    pi_tmp[jj] = 0.0;
  }

  // in_Q gives Q_N = Q_n+1, and Q0[ii] gives Q_0 = Q_n
  // do Runge-Kutta 4 with the 3/8 rule
  for(int mm=0; mm<N; ++mm)
  {
    const double Q_m = Q0[ii] + static_cast<double>(mm) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

    const double Q_mp1 = Q0[ii] + static_cast<double>(mm+1) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

    F_coronary(ii, pi_m, Q_m, dPimdt_k1[ii][mm], K1);

    for(int jj=0; jj<num_odes; ++jj)
      pi_tmp[jj] = pi_m[jj] + fac13 * K1[jj] * h;

    F_coronary(ii, pi_tmp, fac23 * Q_m + fac13 * Q_mp1, dPimdt_k2[ii][mm], K2);

    for(int jj=0; jj<num_odes; ++jj)
      pi_tmp[jj] = pi_m[jj] - fac13*K1[jj] * h + K2[jj] * h;

    F_coronary(ii, pi_tmp, fac13 * Q_m + fac23 * Q_mp1, dPimdt_k3[ii][mm], K3);

    for(int jj=0; jj<num_odes; ++jj)
      pi_tmp[jj] = pi_m[jj] + K1[jj] * h - K2[jj] * h + K3[jj] * h;

    F_coronary(ii, pi_tmp, Q_mp1, dPimdt_k1[ii][mm+1], K4);

    for(int jj=0; jj<num_odes; ++jj)
      pi_m[jj] = pi_m[jj] + fac18 * K1[jj] * h + fac38 * K2[jj] * h + fac38 * K3[jj] * h + fac18 * K4[jj] * h;
  }

  // Make a copy of the ODE solutions.
  // prev_0D_sol will reset initial values Pi0 for future time integration (t=n+1-> t=n+2)
  for(int jj=0; jj<num_odes; ++jj) prev_0D_sol[ii][jj]=pi_m[jj];

  output_P = pi_m[0] + Ra[ii] * in_Q;

  return output_P;
}

double GenBC_Coronary::get_P0( const int &ii ) const
{
  return Q0[ii] * Ra[ii] + Pi0[ii][0];
}

void GenBC_Coronary::reset_initial_sol( const int &ii, const double &in_Q_0,
    const double &in_P_0, const double &curr_time, const bool &is_restart )
{
  Q0[ii] = in_Q_0;

  // For a restart job, use 0D solutions from user-defined time step
  // as initial solutions for the next time integration. Otherwise, use
  // the previous 0D solutions as initial solutions.
  if( is_restart )
    for( int jj=0; jj<num_odes; ++jj ) Pi0[ii][jj] = restart_0D_sol[ii][jj];
  else
    for( int jj=0; jj<num_odes; ++jj ) Pi0[ii][jj] = prev_0D_sol[ii][jj];

  // Precalculate dPimdt values needed for integrating coronary ODEs.
  get_dPim_dt(ii, curr_time, curr_time + N * h);
}

void GenBC_Coronary::F_coronary( const int &ii, const double * const &pi, const double &q,
    const double &dPimdt, double * const &K ) const
{
  // The Coronary LPM consists of two ODEs.

  K[0] = (q-(pi[0]-pi[1])/Ra_micro[ii])/Ca[ii];
  K[1] = ((pi[0]-pi[1])/Ra_micro[ii]-(pi[1]-Pd[ii])/Rv[ii])/Cim[ii]+dPimdt;
}

void GenBC_Coronary::get_dPim_dt( const int &ii, const double &time_start, const double &time_end )
{
  double tend_mod = fmod( time_end, Time_data[ii][num_Pim_data[ii]-1] );
  if ( tend_mod<absTol ) tend_mod =  Time_data[ii][num_Pim_data[ii]-1];

  // Find the interval in Pim that covers current integration time.
  double x1, x2, f1, f2, d1, d2;
  bool is_interval_found = false;
  for(int mm=1; mm<num_Pim_data[ii]; ++mm)
  {
    if ( tend_mod <= Time_data[ii][mm] )
    {
      x1 = Time_data[ii][mm-1];
      f1 = Pim_data[ii][mm-1];
      d1 = der_Pim_data[ii][mm-1];
      x2 = Time_data[ii][mm];
      f2 = Pim_data[ii][mm];
      d2 = der_Pim_data[ii][mm];
      is_interval_found = true;
      break;
    }
  }

  SYS_T::print_fatal_if( !is_interval_found, "Error: get_dPim_dt the interval is not found.\n" );

  std::vector<double> xe_1,xe_2,xe_3;
  xe_1.resize( N+1 );
  xe_2.resize( N );
  xe_3.resize( N );

  // Precalculate dPimdt values at time points evaluated by RK4
  // from time_start to time_end.
  double tmp = time_start;

  const double fac13 = 1.0 / 3.0;
  const double fac23 = 2.0 / 3.0;

  for(int mm=0;mm<N;++mm)
  {
    xe_1[mm] = fmod(tmp,Time_data[ii][num_Pim_data[ii]-1]);
    xe_2[mm] = fmod(tmp+fac13*h,Time_data[ii][num_Pim_data[ii]-1]);
    xe_3[mm] = fmod(tmp+fac23*h,Time_data[ii][num_Pim_data[ii]-1]);
    tmp = tmp + h;
  }

  xe_1[N] = tend_mod;

  GENBC_T::get_cubic_hermite_der( x1, x2, f1,f2, d1, d2, N+1, xe_1, dPimdt_k1[ii] );

  GENBC_T::get_cubic_hermite_der( x1, x2, f1,f2, d1, d2, N,   xe_2, dPimdt_k2[ii] );

  GENBC_T::get_cubic_hermite_der( x1, x2, f1,f2, d1, d2, N,   xe_3, dPimdt_k3[ii] );
}

// EOF
