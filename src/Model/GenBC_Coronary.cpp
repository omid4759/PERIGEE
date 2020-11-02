#include "GenBC_Coronary.hpp"

GenBC_Coronary::GenBC_Coronary( const char * const &lpn_filename, const int &in_N,
    const double &dt3d )
: N( in_N ), h( dt3d/static_cast<double>(N) ),
  absTol( 1.0e-8 ), relTol( 1.0e-5 )
{
  // Now read the lpn files for num_ebc, Ra, Ca, Ra_micro, cim, Rv and Pd
  std::string temp_name( lpn_filename );
  SYS_T::file_check( temp_name ); // make sure the file is on the disk

  std::ifstream reader;
  reader.open( lpn_filename, std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;
  std::string bc_type;

  tstart=0.0;
  tend=N*h;
  // The first non-commented lien should be
  // RCR num_ebc
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

  // Check the file's bc_type matches RCR
  if( bc_type.compare("Coronary") ==0
      || bc_type.compare("CORONARY") == 0
      || bc_type.compare("coronary") == 0 )
  {
    Ra.resize( num_ebc ); Ca.resize( num_ebc );
    Ra_micro.resize( num_ebc ); Cim.resize( num_ebc );
    Rv.resize( num_ebc ); Pd.resize( num_ebc ); alpha_Pim.resize(num_ebc);
    Q0.resize(num_ebc);
    Pi0.resize(num_ebc); num_Pimdata.resize(num_ebc);
    tdata.resize(num_ebc);Pimdata.resize(num_ebc);
    Pimderdata.resize(num_ebc);prev_0D_sol.resize(num_ebc);
    dPimdt_k1.resize(num_ebc);dPimdt_k2.resize(num_ebc);dPimdt_k3.resize(num_ebc);
    for(int ii =0;ii<2;++ii){
      prev_0D_sol[ii].resize(2);
      Pi0[ii].resize(2);
    }

    for(int ii=0;ii<num_ebc;++ii){
      dPimdt_k1[ii].resize(N+1);
      dPimdt_k2[ii].resize(N);
      dPimdt_k3[ii].resize(N);

    }


  }
  else SYS_T::print_fatal("Error: the outflow model in %s does not match GenBC_Coronary.\n", lpn_filename);

  // Read files for each ebc to set the values of Ra, Ca, Ra_micro, Cim,Rv, and Pd
  int counter = 0;
  int data_size=1;
  while( std::getline(reader, sline) )
  {
    if( sline[0] != '#' && !sline.empty() )
    {
      sstrm.str( sline );
      int face_id;
      sstrm >> face_id;

      // Make sure the face_id, the first column in the file are listed
      // from 0 to ebc_id - 1
      if(face_id != counter) SYS_T::print_fatal("Error: GenBC_Coronary the input file %s has wrong format in the face id column (the first column). \n", lpn_filename);

      sstrm >> Ra[ counter ];
      sstrm >> Ca[ counter ];
      sstrm >> Ra_micro[ counter ];
      sstrm >> Cim[ counter ];
      sstrm >> Rv[ counter ];
      sstrm >> Pd[ counter ];
      sstrm >> num_Pimdata[counter];
      sstrm >> alpha_Pim[counter];
      printf("reading coronary outlet=%d Ra=%lf, Ca=%lf, Ramicro=%lf Cim=%lf RV=%lf Pd=%lf numPim=%d \n",counter,Ra[counter],Ca[counter],Ra_micro[counter],Cim[counter],Rv[counter],Pd[counter],num_Pimdata[counter]);
      SYS_T::print_fatal_if(num_Pimdata[counter]<=2 && num_Pimdata[counter]!=0, "Error: num of Pim data needs to be 0 for RCR or >2 for coronary  \n");
      if(num_Pimdata[counter]>0){
       data_size=num_Pimdata[counter];
      }else{
       data_size=1;
      }

      tdata[counter].resize(data_size);
      Pimdata[counter].resize(data_size);
      Pimderdata[counter].resize(data_size);

      sstrm.clear();

      for (int ii =0;ii<num_Pimdata[counter];++ii){
        getline(reader, sline);
        sstrm.str( sline );
        sstrm>>tdata[counter][ii];

        sstrm>>Pimdata[counter][ii];
        Pimdata[counter][ii]=Pimdata[counter][ii]*alpha_Pim[counter];
        sstrm.clear();
        printf("Pim t=%lf Pim*alpha_Pim=%lf \n",tdata[counter][ii],Pimdata[counter][ii]);

      }
      if(num_Pimdata[counter]>0){
        set_phcip(counter);
        get_dPimdt(counter);

      }

      SYS_T::print_fatal_if(tdata[counter][0]>0.0, "Error: Pim data do not start from 0.\n");

      counter += 1;
    }
  }

  if(counter != num_ebc ) SYS_T::print_fatal("Error: GenBC_Coronary the input file %s does not contain complete data for outlet faces. \n", lpn_filename);

  reader.close();

  SYS_T::commPrint( "===> GenBC_Coronary data are read in from %s.\n", lpn_filename );

  // Set a zero initial value. They should be reset based on the initial
  // 3D solutions.
  for(int ii=0; ii<num_ebc; ++ii)
  {

    Q0[ii] = 0.0;
    Pi0[ii][0] = 0.0;
    Pi0[ii][1]=0.0;
    prev_0D_sol[ii][0]=0.0;
    prev_0D_sol[ii][1]=0.0;

    // Make sure C and R are nonzero
    SYS_T::print_fatal_if(Ca[ii]==0.0, "Error: GenBC_Coronary Ca cannot be zero.\n");
    SYS_T::print_fatal_if(Cim[ii]==0.0, "Error: GenBC_Coronary Cim cannot be zero.\n");
    SYS_T::print_fatal_if(Ra[ii]==0.0, "Error: GenBC_Coronary Ra cannot be zero.\n");
    SYS_T::print_fatal_if(Ra_micro[ii]==0.0, "Error: GenBC_Coronary Ra_micro cannot be zero.\n");
    SYS_T::print_fatal_if(Rv[ii]==0.0, "Error: GenBC_Coronary Rv cannot be zero.\n");


  }

}

GenBC_Coronary::~GenBC_Coronary()
{}

void GenBC_Coronary::print_info() const
{
  SYS_T::commPrint( "     Coronary model: N = %d, h = %e, num_ebc = %d \n", N, h, num_ebc );

  for(int ii=0; ii<num_ebc; ++ii)
    SYS_T::commPrint( "     ebcid = %d, Ra = %e, Ca = %e, Ra_micro =%e,Rim =%e,Rv =%e, Pd = %e \n", ii, Ra[ii], Ca[ii],Ra_micro[ii],Cim[ii], Rv[ii], Pd[ii] );
}


double GenBC_Coronary::get_m( const int &ii, const double &in_dot_Q,
   const double &in_Q ) const
{
  double diff = std::abs(in_Q) * relTol;

  if( diff < absTol ) diff = absTol;

  const double left  = get_P(ii, in_dot_Q, in_Q + 0.5 * diff);
  const double right = get_P(ii, in_dot_Q, in_Q - 0.5 * diff);

  return (left - right) / diff;
}


double GenBC_Coronary::get_P( const int &ii, const double &in_dot_Q,
   const double &in_Q ) const
{

  const double fac13 = 1.0 / 3.0;
  const double fac23 = 2.0 / 3.0;
  const double fac18 = 1.0 / 8.0;
  const double fac38 = 3.0 / 8.0;

  //num_Pimdata indicates this is a rcr outlet
  if(num_Pimdata[ii]==0){

    double pi_m = Pi0[ii][0]; // Pi_m

  // in_Q gives Q_N = Q_n+1, and Q0[ii] gives Q_0 = Q_n
  for(int mm=0; mm<N; ++mm)
  {
    const double Q_m = Q0[ii] + static_cast<double>(mm) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

    const double Q_mp1 = Q0[ii] + static_cast<double>(mm+1) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

    const double K1 = F(ii, pi_m, Q_m );

    const double K2 = F(ii, pi_m + fac13 * K1 * h, fac23 * Q_m + fac13 * Q_mp1 );

    const double K3 = F(ii, pi_m - fac13 * K1 * h + K2 * h, fac13 * Q_m + fac23 * Q_mp1);

    const double K4 = F(ii, pi_m + K1 * h - K2 * h + K3 * h, Q_mp1);

    pi_m = pi_m + fac18 * K1 * h + fac38 * K2 * h + fac38 * K3 * h + fac18 * K4 * h;
  }
  prev_0D_sol[ii][0]=pi_m;


  return pi_m + Ra[ii] * in_Q ;
  //return pi_m + Ra[ii] * in_Q + Pd[ii];

  }


  double pi_m[2];

  //initial Pressures at Ca and Cim
  pi_m[0]= Pi0[ii][0];
  pi_m[1]= Pi0[ii][1];

  // in_Q gives Q_N = Q_n+1, and Q0[ii] gives Q_0 = Q_n



  double K1[2],K2[2],K3[2],K4[2];
  double pi_tmp[2];


  for(int mm=0; mm<N; ++mm)
  {
    const double Q_m = Q0[ii] + static_cast<double>(mm) * ( in_Q - Q0[ii] ) / static_cast<double>(N);

    const double Q_mp1 = Q0[ii] + static_cast<double>(mm+1) * ( in_Q - Q0[ii] ) / static_cast<double>(N);


   F(ii, pi_m, Q_m, dPimdt_k1[ii][mm],K1 );


    for(int j =0; j<2; ++j){
    pi_tmp[j]=pi_m[j]+fac13*K1[j]*h;
    }

    F(ii, pi_tmp , fac23 * Q_m + fac13 * Q_mp1,dPimdt_k2[ii][mm],K2);


    for(int j =0; j<2; ++j){
    pi_tmp[j]=pi_m[j]-fac13*K1[j]*h+K2[j]*h;

    }


    F(ii, pi_tmp , fac13 * Q_m + fac23 * Q_mp1,dPimdt_k3[ii][mm],K3 );


    for(int j =0; j<2; ++j){
    pi_tmp[j]=pi_m[j] + K1[j] * h - K2[j] * h + K3[j] * h;
    }


    F(ii, pi_tmp, Q_mp1, dPimdt_k1[ii][mm+1],K4);



    for(int j =0; j<2; ++j){
      pi_m[j] = pi_m[j] + fac18 * K1[j] * h + fac38 * K2[j] * h + fac38 * K3[j] * h + fac18 * K4[j] * h;
    }
  }



  prev_0D_sol[ii][0]=pi_m[0];
  prev_0D_sol[ii][1]=pi_m[1];


  return pi_m[0] + Ra[ii] * in_Q;

  //return pi_m[0] + Ra[ii] * in_Q + Pd[ii];

}


void GenBC_Coronary::reset_initial_sol( const int &ii, const double &in_Q_0,
    const double &in_P_0, const double &curr_time )
{
  Q0[ii]  = in_Q_0;

  //when reset_initial_sol is called, use the last converged 0D solition as initial solutions.
  Pi0[ii][0] = prev_0D_sol[ii][0];
  //Pi0[ii][0] = in_P_0-in_Q_0*Ra[ii];

  Pi0[ii][1] = prev_0D_sol[ii][1];


  //Pi0[ii][0] = in_P_0 - in_Q_0 * Ra[ii] - Pd[ii];


  if(ii==0){
  tstart=curr_time;
  tend=curr_time+N*h;
   }

  if(num_Pimdata[ii]>2){
     get_dPimdt(ii);
    }


}

double GenBC_Coronary:: F(const int &ii, const double *pi, const double &q, const double &dPimdt, double * K)const
{
      //return -1.0 * pi / (Rd[ii] * C[ii]) + q / C[ii];
      K[0]=(q-(pi[0]-pi[1])/Ra_micro[ii])/Ca[ii];
      K[1]=((pi[0]-pi[1])/Ra_micro[ii]-(pi[1]-Pd[ii])/(Rv[ii]))/Cim[ii]+dPimdt;

    //  K[0]=(q-(pi[0]-Pd[ii])/Ra_micro[ii])/Ca[ii];
    //  K[1]=((pi[0]-pi[1])/Ra_micro[ii]-(pi[1]-Pd[ii])/(Rv[ii]))/Cim[ii];

     //when pi is the pressure difference across the capacitor
    //  K[0]=(q-(pi[0]-pi[1])/Ra_micro[ii])/Ca[ii];
    //  K[1]=((pi[0]-pi[1])/Ra_micro[ii]-(pi[1])/(Rv[ii]))/Cim[ii]+dPimdt;


    //  K[0]=-1.0 * pi[0] / (Ra_micro[ii] * Ca[ii]) + q / Ca[ii];
      return 0.0;
}


double GenBC_Coronary:: F(const int &ii, const double &pi, const double &q) const
{
    // return -1.0 * pi / (Ra_micro[ii] * Ca[ii]) + q / Ca[ii];

     return (q-(pi-Pd[ii])/Ra_micro[ii])/Ca[ii];
}


void  GenBC_Coronary:: set_phcip(const int &ii)
{

   spline_pchip_set ( num_Pimdata[ii], tdata[ii], Pimdata[ii], Pimderdata[ii] );

}


void GenBC_Coronary:: spline_pchip_set (int n, std::vector<double> &x, std::vector<double> &f, std::vector<double> &d)
//****************************************************************************80
//
//  Purpose:
//
//    SPLINE_PCHIP_SET sets derivatives for a piecewise cubic Hermite interpolant.
//
//  Discussion:
//
//    This routine computes what would normally be called a Hermite
//    interpolant.  However, the user is only required to supply function
//    values, not derivative values as well.  This routine computes
//    "suitable" derivative values, so that the resulting Hermite interpolant
//    has desirable shape and monotonicity properties.
//
//    The interpolant will have an extremum at each point where
//    monotonicity switches direction.
//
//    The resulting piecewise cubic Hermite function may be evaluated
//    by SPLINE_PCHIP_VAL..
//
//    This routine was originally called "PCHIM".
//
//    An "abs" was corrected to a "fabs" on the report of Thomas Beutlich,
//    10 October 2012.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2005
//
//  Author:
//
//    FORTRAN77 original version by Fred Fritsch, Lawrence Livermore National Laboratory.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//    Fred Fritsch, Judy Butland,
//    A Method for Constructing Local Monotone Piecewise
//    Cubic Interpolants,
//    SIAM Journal on Scientific and Statistical Computing,
//    Volume 5, Number 2, 1984, pages 300-304.
//
//  Parameters:
//
//    Input, int N, the number of data points.  N must be at least 2.
//
//    Input, double X[N], the strictly increasing independent
//    variable values.
//
//    Input, double F[N], dependent variable values to be interpolated.  This
//    routine is designed for monotonic data, but it will work for any F-array.
//    It will force extrema at points where monotonicity switches direction.
//
//    Output, double D[N], the derivative values at the
//    data points.  If the data are monotonic, these values will determine
//    a monotone cubic Hermite function.
//
{
  double del1;
  double del2;
  double dmax;
  double dmin;
  double drat1;
  double drat2;
  double dsave;
  double h1;
  double h2;
  double hsum;
  double hsumt3;
  int i;
  int ierr;
  int nless1;
  double temp;
  double w1;
  double w2;
//
//  Check the arguments.
//

  SYS_T::print_fatal_if(n<2, "Error: GenBC_Coronary SPLINE_PCHIP_SET: Number of evaluation points is less than 1 \n");

  for ( i = 1; i < n; i++ )
  {
      SYS_T::print_fatal_if(x[i] <= x[i-1], "Error: GenBC_Coronary SPLINE_PCHIP_SET: X array not strictly increasing. \n");

  }

  ierr = 0;
  nless1 = n - 1;
  h1 = x[1] - x[0];
  del1 = ( f[1] - f[0] ) / h1;
  dsave = del1;
//
//  Special case N=2, use linear interpolation.
//
  if ( n == 2 )
  {
    d[0] = del1;
    d[n-1] = del1;
    return;
  }
//
//  Normal case, 3 <= N.
//
  h2 = x[2] - x[1];
  del2 = ( f[2] - f[1] ) / h2;
//
//  Set D(1) via non-centered three point formula, adjusted to be
//  shape preserving.
//
  hsum = h1 + h2;
  w1 = ( h1 + hsum ) / hsum;
  w2 = -h1 / hsum;
  d[0] = w1 * del1 + w2 * del2;

  if ( pchst ( d[0], del1 ) <= 0.0 )
  {
    d[0] = 0.0;
  }
//
//  Need do this check only if monotonicity switches.
//
  else if ( pchst ( del1, del2 ) < 0.0 )
  {
     dmax = 3.0 * del1;

     if ( fabs ( dmax ) < fabs ( d[0] ) )
     {
       d[0] = dmax;
     }

  }
//
//  Loop through interior points.
//
  for ( i = 2; i <= nless1; i++ )
  {
    if ( 2 < i )
    {
      h1 = h2;
      h2 = x[i] - x[i-1];
      hsum = h1 + h2;
      del1 = del2;
      del2 = ( f[i] - f[i-1] ) / h2;
    }
//
//  Set D(I)=0 unless data are strictly monotonic.
//
    d[i-1] = 0.0;

    temp = pchst ( del1, del2 );

    if ( temp < 0.0 )
    {
      ierr = ierr + 1;
      dsave = del2;
    }
//
//  Count number of changes in direction of monotonicity.
//
    else if ( temp == 0.0 )
    {
      if ( del2 != 0.0 )
      {
        if ( pchst ( dsave, del2 ) < 0.0 )
        {
          ierr = ierr + 1;
        }
        dsave = del2;
      }
    }
//
//  Use Brodlie modification of Butland formula.
//
    else
    {
      hsumt3 = 3.0 * hsum;
      w1 = ( hsum + h1 ) / hsumt3;
      w2 = ( hsum + h2 ) / hsumt3;
      dmax = r8_max ( fabs ( del1 ), fabs ( del2 ) );
      dmin = r8_min ( fabs ( del1 ), fabs ( del2 ) );
      drat1 = del1 / dmax;
      drat2 = del2 / dmax;
      d[i-1] = dmin / ( w1 * drat1 + w2 * drat2 );
    }
  }
//
//  Set D(N) via non-centered three point formula, adjusted to be
//  shape preserving.
//
  w1 = -h2 / hsum;
  w2 = ( h2 + hsum ) / hsum;
  d[n-1] = w1 * del1 + w2 * del2;

  if ( pchst ( d[n-1], del2 ) <= 0.0 )
  {
    d[n-1] = 0.0;
  }
  else if ( pchst ( del1, del2 ) < 0.0 )
  {
//
//  Need do this check only if monotonicity switches.
//
    dmax = 3.0 * del2;

    if ( fabs ( dmax ) < fabs ( d[n-1] ) )
    {
      d[n-1] = dmax;
    }

  }
  return;


}

double GenBC_Coronary:: pchst ( double arg1, double arg2 )

//****************************************************************************80
//
//  Purpose:
//
//    PCHST: PCHIP sign-testing routine.
//
//  Discussion:
//
//    This routine essentially computes the sign of ARG1 * ARG2.
//
//    The object is to do this without multiplying ARG1 * ARG2, to avoid
//    possible over/underflow problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 August 2005
//
//  Author:
//
//    Original FORTRAN77 version by Fred Fritsch, Lawrence Livermore National Laboratory.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Fred Fritsch, Ralph Carlson,
//    Monotone Piecewise Cubic Interpolation,
//    SIAM Journal on Numerical Analysis,
//    Volume 17, Number 2, April 1980, pages 238-246.
//
//  Parameters:
//
//    Input, double ARG1, ARG2, two values to check.
//
//    Output, double PCHST,
//    -1.0, if ARG1 and ARG2 are of opposite sign.
//     0.0, if either argument is zero.
//    +1.0, if ARG1 and ARG2 are of the same sign.
//
{
  double value;

  if ( arg1 == 0.0 )
  {
    value = 0.0;
  }
  else if ( arg1 < 0.0 )
  {
    if ( arg2 < 0.0 )
    {
      value = 1.0;
    }
    else if ( arg2 == 0.0 )
    {
      value = 0.0;
    }
    else if ( 0.0 < arg2 )
    {
      value = -1.0;
    }
  }
  else if ( 0.0 < arg1 )
  {
    if ( arg2 < 0.0 )
    {
      value = -1.0;
    }
    else if ( arg2 == 0.0 )
    {
      value = 0.0;
    }
    else if ( 0.0 < arg2 )
    {
      value = 1.0;
    }
  }

  return value;
}

double GenBC_Coronary::r8_max ( double x, double y )const

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  if ( y < x )
  {
    return x;
  }
  else
  {
    return y;
  }
}
//****************************************************************************80

double GenBC_Coronary::r8_min ( double x, double y )const

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  if ( y < x )
  {
    return y;
  }
  else
  {
    return x;
  }
}

void GenBC_Coronary:: get_dPimdt(const int &ii)
{

  double x1,x2,f1,f2,d1,d2;

  const double fac13 = 1.0 / 3.0;
  const double fac23 = 2.0 / 3.0;
//  const PetscMPIInt rank = SYS_T::get_MPI_rank();

  double tend_mod=fmod(tend,tdata[ii][num_Pimdata[ii]-1]);
  if (fabs(tend_mod)<1e-7){
    tend_mod=tdata[ii][num_Pimdata[ii]-1];
  }

     for(int mm=1; mm<num_Pimdata[ii];++mm){
       if (tend_mod<=tdata[ii][mm]){

        x1=tdata[ii][mm-1];
        f1=Pimdata[ii][mm-1];
        d1=Pimderdata[ii][mm-1];
        x2=tdata[ii][mm];
        f2=Pimdata[ii][mm];
        d2=Pimderdata[ii][mm];

        break;
       }
    }



  std::vector<double> xe_1,xe_2,xe_3;
  xe_1.resize( N+1 );
  xe_2.resize(N);
  xe_3.resize(N);

   double tmp=tstart;

/*
  if(rank == 0 ){
   printf("check time mod(time) tstart=%lf tend=%lf N=%d x1=%lf f1=%lf d1=%lf x2=%lf f2=%lf d2=%lf \n",tstart,tend,N,x1,f1,d1,x2,f2,d2);
   }
   */
  for (int mm=0;mm<N;++mm){



    xe_1[mm]=fmod(tmp,tdata[ii][num_Pimdata[ii]-1]);
    xe_2[mm]=fmod(tmp+fac13*h,tdata[ii][num_Pimdata[ii]-1]);
    xe_3[mm]=fmod(tmp+fac23*h,tdata[ii][num_Pimdata[ii]-1]);


    tmp=tmp+h;

  }

  xe_1[N]=tend_mod;


  cubic_hermite_derivative ( x1, x2, f1,f2, d1, d2, N+1, xe_1, dPimdt_k1[ii]);

  cubic_hermite_derivative ( x1, x2, f1,f2, d1, d2, N, xe_2, dPimdt_k2[ii]);

  cubic_hermite_derivative ( x1, x2, f1,f2, d1, d2, N, xe_3, dPimdt_k3[ii]);
/*
  if(rank == 0){
   printf("check  mod(time) dPimdt \n");
  for (int kk=0;kk<N;kk++){
     printf("%lf %lf \n",xe_1[kk],dPimdt_k1[kk]);
     printf("%lf %lf \n",xe_2[kk],dPimdt_k2[kk]);
     printf("%lf %lf \n",xe_3[kk],dPimdt_k3[kk]);
  }

    printf("%lf %lf \n",xe_1[N],dPimdt_k1[N]);
    }
*/


}


void GenBC_Coronary:: cubic_hermite_derivative ( double x1, double x2, double f1, double f2, double d1, double d2,
  int ne,std::vector<double> &xe, std::vector<double> &fe )

//****************************************************************************80
// hermite cubic interpolation

{

  double t;
  double h;
  int i;
  double c2;
  double c3;


  SYS_T::print_fatal_if(ne<1, "Error: GenBC_Coronary cubic_hermite_derivative: Number of evaluation points is less than 1 \n");

  h = x2 - x1;

  SYS_T::print_fatal_if(h==0.0, "Error: GenBC_Coronary cubic_hermite_derivative: The interval [X1,X2] is of zero length \n");
//
//  Initialize.
//



   c3=6.0*f1/h+3.0*d1-6.0*f2/h+3.0*d2;

   c2=-6.0*f1/h-4.0*d1+6.0*f2/h-2.0*d2;



//  Evaluation loop.
//



  for ( i = 0; i < ne; i++ )
  {
    t = (xe[i] - x1)/h;

    fe[i] = t * ( c2 + t * c3 )+d1 ;

  }



}

// EOF
