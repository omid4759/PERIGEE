#ifndef GENBC_Coronary_HPP
#define GENBC_Coronary_HPP
// ==================================================================
// GenBC_Coronary.hpp
// This files defines the General boundary condition in cardiovascular
// flow simulations.
// wgyang 2020/10
// ==================================================================
#include "IGenBC.hpp"

class GenBC_Coronary : public IGenBC
{
  public:
    GenBC_Coronary( const char * const &lpn_filename, const int &in_N,
       const double &dt3d );

    virtual ~GenBC_Coronary();

    virtual void print_info() const;

    virtual int get_num_ebc() const {return num_ebc;}

    virtual double get_m( const int &ii, const double &in_dot_Q,
       const double &in_Q ) const ;

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
    //  return Q0[ii] * Ra[ii] + Pi0[ii][0] +Pd[ii] ;
      return Q0[ii] * Ra[ii] + Pi0[ii][0]  ;

    }

    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
        const double &in_P_0, const double & curr_time );

    virtual void get_m( double * const &in_dot_Q,
       double * const &in_Q, double * const &m ) const ;

    virtual void get_n( double * const &in_dot_Q,
       double * const &in_Q, double * const &n ) const
    {
      for(int ii=0;ii<num_ebc;++ii){
       n[ii]=0.0;
      }
    }

    // Obtain P in order to the define the outlet traction for the ii-th
    // outlet surface

    virtual void get_P( double * const &in_dot_Q,
       double * const &in_Q, double * const &P ) const;


    virtual void get_P0( double * const &Pn ) const
    {
      for (int ii=0;ii<num_ebc;++ii){
       Pn[ii]= Q0[ii] * Ra[ii] + Pi0[ii][0]  ;
      }
    }

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

    // Vectors storing the Ra, Ca, Ra_micro, Cim, Rv and Pd values on each outlet faces
    // the length of the following 4 vectors are num_ebc
    std::vector<double> Ra, Ca, Ra_micro,Cim, Rv, Pd,alpha_Pim; // R-C-R-C-R coronary parameters
    std::vector<std::vector<double> > tdata,Pimdata,Pimderdata, Pi0;//time and pressure for Pim
    std::vector<std::vector<double>> dPimdt_k1,dPimdt_k2,dPimdt_k3;
    mutable std::vector<std::vector<double> > prev_0D_sol;


    std::vector<int> num_Pimdata;

    // Vectors storing the Q0 and Pi0 on each outlet faces
    std::vector<double> Q0;


    double F(const int &ii, const double *pi, const double &q, const double &dPimdt, double * K)const;//for coronary outlets
    double F(const int &ii, const double &pi, const double &q) const; //for rcr outlets
    void get_dPimdt(const int &ii ) ;
    void  cubic_hermite_derivative ( double x1, double x2, double f1, double f2, double d1, double d2,
      int ne, std::vector<double> &xe, std::vector<double> &fe );
    void set_phcip(const int &ii);
    void spline_pchip_set (int n, std::vector<double> &x, std::vector<double> &f, std::vector<double> &d);
    double pchst ( double arg1, double arg2 );
    double r8_max ( double x, double y )const ;
    double r8_min ( double x, double y )const;

};

#endif
