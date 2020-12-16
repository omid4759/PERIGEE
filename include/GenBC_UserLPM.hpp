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
    virtual double get_P0( const int &ii ) const
    {
      return P0_Neumann[ii];
    }

    virtual void reset_initial_sol( double * const &in_Q_0,
       double * const &in_P_0, const double & curr_time );

//  for GenBC with Dirichlet faces
    virtual void reset_initial_sol( double * const &in_Q_0_Neumann,
       double * const &in_P_0_Neumann, double * const &in_Q_0_Dirichlet,
       double * const &in_P_0_Dirichlet, const double &curr_time ) ;

    virtual void get_m( double * const &in_dot_Q,
       double * const &in_Q, double * const &in_P, double * const &m ) const ;

    virtual double get_curr_P(const int &ii )const;

    virtual double get_curr_Q(const int &ii)const;

    virtual double get_curr_m(const int &ii)const;

    virtual double get_curr_n(const int &ii)const;


    virtual void set_curr_P(const int & ii, const double & Pi);

    virtual void set_curr_Q(const int & ii, const double & Qi);

    virtual void set_curr_m(const int & ii,const double & mi);

    virtual void set_curr_n(const int &ii, const double &ni);

    virtual void get_P_Q(double * const &in_dot_Q,
       double * const &in_Q, double * const &in_P, double * const &P_Neumann, double * const &Q_Dirichlet, const bool & output_alldata_flag)const;

    virtual void get_Q0(double * const &Qn) const;
    double get_Q0( const int &ii ) const;
    virtual int get_num_Dirichlet_faces()const {return num_Dirichlet_faces;}



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

    std::vector<int> Neumann_Surface_To_LPM;
    std::vector<int> Dirichlet_Surface_To_LPM;

    int num_Dirichlet_faces;

    mutable std::vector<double> prev_0D_sol,Xi0;
    std::vector<std::string> face_type;


    // Vectors storing the Q0 and Pi0 on each Neumann and Dirichlet faces respectively.
    std::vector<double> Q0_Neumann,P0_Neumann,Q0_Dirichlet,P0_Dirichlet;

    // Vectors storing the current P and Q on each Neumann and Dirichlet faces respectively after complete of ODE integration.
    std::vector<double> curr_P,curr_Q,curr_m,curr_n;

};

#endif
