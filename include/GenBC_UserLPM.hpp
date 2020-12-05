#ifndef GENBC_UserLPM_HPP
#define GENBC_UserLPM_HPP
// ==================================================================
// GenBC_Coronary.hpp
// This files defines the General boundary condition in cardiovascular
// flow simulations.
// wgyang 2020/10
// ==================================================================
#include "IGenBC.hpp"

#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>



class GenBC_UserLPM : public IGenBC
{
  public:
    GenBC_UserLPM( const char * const &lpn_filename, const int &in_N,
       const double &dt3d );

    virtual ~GenBC_UserLPM();

    virtual void print_info() const;

    virtual int get_num_ebc() const {return num_ebc;}

    virtual void get_m( double * const &in_dot_Q,
       double * const &in_Q, double * const &m ) const ;

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


    virtual void get_P0(double * const &Pn) const;


    virtual void reset_initial_sol( double * const &in_Q_0,
       double * const &in_P_0, const double & curr_time );

  private:
    const int N;

    const double h; // delta t = Nh

    // Parameters used to define difference quotient for get_m.
    const double absTol, relTol;

    int num_ebc;

    double tstart;
    double tend;
    int num_LPM_unknowns;

    const char * myfifo1 = "/tmp/myfifo1";
    const char * myfifo2 = "/tmp/myfifo2";

    std::vector<int> Surface_To_LPM;

    mutable std::vector<double> prev_0D_sol,Xi0;
    std::vector<std::string> face_type;


    // Vectors storing the Q0 and Pi0 on each outlet faces
    std::vector<double> Q0,P0;


    double F(double * const &pi, double * const &q, double * const &K)const;
};

#endif
