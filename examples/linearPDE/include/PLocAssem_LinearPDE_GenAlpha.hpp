#ifndef PLOCASSEM_LINEARPDE_GENALPHA_HPP
#define PLOCASSEM_LINEARPDE_GENALPHA_HPP
// ============================================================================
// PLocAssem_LinearPDE_GenAlpha.hpp
//
// Date: Oct 6 2023
// ============================================================================
#include "IPLocAssem_Linear.hpp"
#include "TimeMethod_GenAlpha.hpp"

class PLocAssem_LinearPDE_GenAlpha : public IPLocAssem_Linear
{
  public:
    PLocAssem_LinearPDE_GenAlpha( 
        const double &in_rho, const double &in_cap, const double &in_kappa,
        const TimeMethod_GenAlpha * const &tm_gAlpha,
        const int &in_nlocbas, const int &in_snlocbas,
        const int &in_num_ebc_fun, const int &in_dof,
        const int &in_dof_mat, const int &elemtype = 501 );

    virtual ~PLocAssem_LinearPDE_GenAlpha();

    virtual int get_dof() const {return dof;}

    virtual int get_dof_mat() const {return dof_mat;}

    virtual void Zero_Mass_Stiffness_Load()
    {
      for(int ii=0; ii<vec_size; ++ii) Load[ii] = 0.0;
      for(int ii=0; ii<vec_size*vec_size; ++ii) 
      {
        Mass[ii] = 0.0;
        Stiffness[ii] = 0.0;
      }
    }

    virtual void Zero_Load()
    {
      for(int ii=0; ii<vec_size; ++ii) Load[ii] = 0.0;
    }

    virtual void Zero_sur_Mass_Load()
    {
      for(int ii=0; ii<sur_size; ++ii) sur_Load[ii] = 0.0;
      for(int ii=0; ii<sur_size*sur_size; ++ii) sur_Mass[ii] = 0.0;
    }

    virtual void Assem_Estimate()
    {
      for(int ii=0; ii<vec_size*vec_size; ++ii) 
      {
        Mass[ii] = 1.0
        Stiffness[ii] = 1.0;
      }
    }

    virtual void Assem_Load(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Stiffness(
        const double &time, const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Mass(
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

    virtual void Assem_Load_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad );

  private:
    // Private data
    const double rho, cap, kappa;
    const double alpha_f, alpha_m, gamma;
    
    const int num_ebc_fun;
    const int dof, dof_mat;

    int nLocBas, snLocBas, vec_size, sur_size;

    void print_info() const;

    Vector_3 get_f( const Vector_3 &pt, const double &tt ) const
    {
      return Vector_3(0.0, 0.0, 0.0); 
    }

    typedef double ( PLocAssem_LinearPDE_GenAlpha::*locassem_transport_funs )( const Vector_3 &pt, const double &t ) const;

    locassem_transport_funs * flist;

    double get_ebc_fun( const int &ebc_id, const Vector_3 &pt, const double &tt ) const
    {
      return ((*this).*(flist[ebc_id]))(pt, tt);
    }

    double get_g_0( const Vector_3 &pt, const double &time ) const
    {
      return 0.0;
    }

    double get_g_1( const Vector_3 &pt, const double &time ) const
    {
      return 0.0;
    }
};

#endif
