#ifndef GENBC_RCR_HPP
#define GENBC_RCR_HPP
// ==================================================================
// GenBC_RCR.hpp
// This files defines the General boundary condition in cardiovascular
// flow simulations.
//
// Pi := P - Rp Q is the pressure over the capacitor
//
// Author: Ju Liu
// Date: Aug. 19 2019
// ==================================================================
#include "IGenBC.hpp"

class GenBC_RCR : public IGenBC
{
  public:
    GenBC_RCR( const char * const &lpn_filename, const int &in_N,
       const double &dt3d );

    virtual ~GenBC_RCR();

    virtual void print_info() const;

    virtual int get_num_ebc() const {return num_ebc;}

    virtual double get_m( const int &ii, const double &in_dot_Q,
       const double &in_Q ) const;

    virtual double get_n( const int &ii, const double &in_dot_Q,
       const double &in_Q ) const
    {
      return 0.0;
    }

    // Obtain P in order to the define the outlet traction for the ii-th
    // outlet surface
    virtual double get_P( const int &ii, const double &in_dot_Q,
       const double &in_Q ) const;

    virtual double get_P0( const int &ii ) const
    {
      return Q0[ii] * Rp[ii] + Pi0[ii] + Pd[ii];
    }

    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
        const double &in_P_0, const double &curr_time );


    virtual void get_m( double * const &in_dot_Q,
       double *const &in_Q, double * const &m ) const;

    virtual void get_n( double * const &in_dot_Q,
       double * const &in_Q, double * const &n ) const
    {
      for (int ii=0;ii<num_ebc;++ii){
       n[ii]=0.0;
      }
    }

    // Obtain P in order to the define the outlet traction for the ii-th
    // outlet surface
    virtual void get_P( double * const &in_dot_Q,
      double * const &in_Q, double * const &P ) const;

    virtual void get_P0( double * const  & Pn ) const
    {

      for (int ii =0; ii<num_ebc;++ii){
       Pn[ii]= Q0[ii] * Rp[ii] + Pi0[ii] + Pd[ii];
      }

    }

    virtual void reset_initial_sol( double * const &in_Q_0,
       double * const &in_P_0, const double &curr_time );




  private:
    const int N;

    const double h; // delta t = Nh

    // Parameters used to define difference quotient for get_m.
    const double absTol, relTol;

    int num_ebc;

    // Vectors storing the Rp, C, and Rd values on each outlet faces
    // the length of the following 4 vectors are num_ebc
    std::vector<double> Rp, C, Rd, Pd; // R-C-R parameters

    // Vectors storing the Q0 and Pi0 on each outlet faces
    std::vector<double> Q0, Pi0;

    double F(const int &ii, const double &pi, const double &q) const
    {
      return -1.0 * pi / (Rd[ii] * C[ii]) + q / C[ii];
    }



};

#endif
