#include "Post_error_ns.hpp"

double POST_T_NS::exact_pres( const double &x, const double &y, const double &z,
          const double &t )
{
  return k0 * z + std::real( k1 * z * exp(i1*omega*t) ); 
}


void POST_T_NS::exact_velo( const double &x, const double &y, const double &z,
    const double &t, double &val_x, double &val_y, double &val_z )
{
  const double r     = sqrt(x*x+y*y);        // radial coord
  const auto   xi    = Lambda * r / R_pipe;
  const auto bes_top = sp_bessel::besselJ(0, xi);
  const auto bes_bot = sp_bessel::besselJ(0, Lambda);

  val_x = 0.0;
  val_y = 0.0;
  val_z = k0*(x*x + y*y - R_pipe*R_pipe) / (4.0*mu) + std::real( coef1 * exp(i1*omega*t) * (1.0 - bes_top / bes_bot ) );
}


void POST_T_NS::exact_grad_pres( const double &x, const double &y, const double &z,
    const double &t, double &val_x, double &val_y, double &val_z )
{
  val_x = 0.0;
  val_y = 0.0;
  val_z = k0 + std::real( k1 * exp(i1*omega*t) );
}


void POST_T_NS::exact_grad_velo( const double &x, const double &y, const double &z,
    const double &t, Matrix_3x3 &grad_velo )
{
  const double r     = sqrt(x*x+y*y);        // radial coord
  const auto   xi    = Lambda * r / R_pipe;
  const auto bes_top = sp_bessel::besselJ(1, xi);
  const auto bes_bot = sp_bessel::besselJ(0, Lambda);

  grad_velo(0,0) = 0.0;
  grad_velo(0,1) = 0.0;
  grad_velo(0,2) = 0.0;

  grad_velo(1,0) = 0.0;
  grad_velo(1,1) = 0.0;
  grad_velo(1,2) = 0.0;

  grad_velo(2,0) = k0 * x / (2.0*mu) + std::real( coef1 * exp(i1*omega*t) * bes_top * i1_1d5 * Omega * x / (bes_bot * r * R_pipe) );
  grad_velo(2,1) = k0 * y / (2.0*mu) + std::real( coef1 * exp(i1*omega*t) * bes_top * i1_1d5 * Omega * y / (bes_bot * r * R_pipe) );
  grad_velo(2,2) = 0.0;
}


void POST_T_NS::exact_wss( const double &x, const double &y, const double &z,
      const double &t, double &val_x, double &val_y, double &val_z )
{
  const auto bes_top = sp_bessel::besselJ(1, Lambda);
  const auto bes_bot = sp_bessel::besselJ(0, Lambda);

  val_x = 0.0;
  val_y = 0.0;
  val_z = k0 * R_pipe / 2.0 - std::real( k1 * R_pipe * i1_0d5 * exp(i1*omega*t) * bes_top / (Omega * bes_bot) ); 
}


double POST_T_NS::get_pres_l2_error( const double * const &sol,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const IQuadPts * const &quad,
    double * const &R,
    const double &t )
{
  const int nqp = quad -> get_num_quadPts();
  const int nLocBas = element -> get_nLocBas();
  double error = 0.0;

  for(int qua=0; qua<nqp; ++qua)
  {
    double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    double sol_qua = 0.0;
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    element -> get_R(qua, R);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x  += ectrlPts_x[ii] * R[ii];
      coor_y  += ectrlPts_y[ii] * R[ii];
      coor_z  += ectrlPts_z[ii] * R[ii];

      sol_qua += sol[ii] * R[ii];
    }

    double exa_qua = exact_pres( coor_x, coor_y, coor_z, t );
    error += (sol_qua - exa_qua) * (sol_qua - exa_qua) * gwts;
  }

  return error;
}


double POST_T_NS::get_pres_h1_error( const double * const &sol,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const IQuadPts * const &quad,
    double * const &R,
    double * const &Rx, double * const &Ry, double * const &Rz,
    const double &t )
{
  const int nqp = quad -> get_num_quadPts();
  const int nLocBas = element -> get_nLocBas();
  double error = 0.0;
  double exa_x, exa_y, exa_z;

  for(int qua=0; qua<nqp; ++qua)
  {
    double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    double sol_x  = 0.0, sol_y  = 0.0, sol_z  = 0.0;
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    element -> get_R_gradR(qua, R, Rx, Ry, Rz);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];

      sol_x  += sol[ii] * Rx[ii];
      sol_y  += sol[ii] * Ry[ii];
      sol_z  += sol[ii] * Rz[ii];
    }

    exact_grad_pres( coor_x, coor_y, coor_z, t, exa_x, exa_y, exa_z );

    error += (sol_x - exa_x) * (sol_x - exa_x) * gwts;
    error += (sol_y - exa_y) * (sol_y - exa_y) * gwts;
    error += (sol_z - exa_z) * (sol_z - exa_z) * gwts;
  }

  return error;
}


double POST_T_NS::get_velo_l2_error( const double * const &solu,
    const double * const &solv, const double * const &solw,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const IQuadPts * const &quad,
    double * const &R,
    const double &t )
{
  const int nqp = quad -> get_num_quadPts();
  const int nLocBas = element -> get_nLocBas();
  double error = 0.0;
  double exa_u, exa_v, exa_w;

  for(int qua=0; qua<nqp; ++qua)
  {
    double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    double sol_u  = 0.0, sol_v  = 0.0, sol_w  = 0.0;
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    element -> get_R(qua, R);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];

      sol_u  += solu[ii] * R[ii];
      sol_v  += solv[ii] * R[ii];
      sol_w  += solw[ii] * R[ii];
    }

    exact_velo( coor_x, coor_y, coor_z, t, exa_u, exa_v, exa_w );

    error += (sol_u - exa_u) * (sol_u - exa_u) * gwts;
    error += (sol_v - exa_v) * (sol_v - exa_v) * gwts;
    error += (sol_w - exa_w) * (sol_w - exa_w) * gwts;
  }

  return error;
}


double POST_T_NS::get_velo_h1_error( const double * const &solu,
    const double * const &solv, const double * const &solw,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const IQuadPts * const &quad,
    double * const &R,
    double * const &Rx, double * const &Ry, double * const &Rz,
    const double &t )
{
  const int nqp = quad -> get_num_quadPts();
  const int nLocBas = element -> get_nLocBas();
  double error = 0.0;

  Matrix_3x3 exa_grad_velo; exa_grad_velo.gen_zero();
  Matrix_3x3 sol_grad_velo; sol_grad_velo.gen_zero();

  for(int qua=0; qua<nqp; ++qua)
  {
    double gwts = element->get_detJac(qua) * quad->get_qw(qua);
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;
    sol_grad_velo.gen_zero();

    element -> get_R_gradR(qua, R, Rx, Ry, Rz);
    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];

      sol_grad_velo(0,0) += solu[ii] * Rx[ii];
      sol_grad_velo(0,1) += solu[ii] * Ry[ii];
      sol_grad_velo(0,2) += solu[ii] * Rz[ii];

      sol_grad_velo(1,0) += solv[ii] * Rx[ii];
      sol_grad_velo(1,1) += solv[ii] * Ry[ii];
      sol_grad_velo(1,2) += solv[ii] * Rz[ii];

      sol_grad_velo(2,0) += solw[ii] * Rx[ii];
      sol_grad_velo(2,1) += solw[ii] * Ry[ii];
      sol_grad_velo(2,2) += solw[ii] * Rz[ii];
    }

    exact_grad_velo( coor_x, coor_y, coor_z, t, exa_grad_velo );

    for(int ii=0; ii<9; ++ii) error += ( exa_grad_velo(ii) - sol_grad_velo(ii) ) * ( exa_grad_velo(ii) - sol_grad_velo(ii) ) * gwts;
  }

  return error;
}


double POST_T_NS::get_wss_l2_error( const double * const &solu,
      const double * const &solv, const double * const &solw,
      const FEAElement * const &element_s,
      const double * const &sctrlPts_x,
      const double * const &sctrlPts_y,
      const double * const &sctrlPts_z,
      const IQuadPts * const &quad_s,
      double * const &R,
      const double &t )
{
  const int nqp = quad_s -> get_num_quadPts();
  const int snLocBas = element_s -> get_nLocBas();
  double error = 0.0;
  double exa_u, exa_v, exa_w;

  for(int qua=0; qua<nqp; ++qua)
  {
    double gwts = element_s->get_detJac(qua) * quad_s->get_qw(qua);
    double sol_u  = 0.0, sol_v  = 0.0, sol_w  = 0.0;
    double coor_x = 0.0, coor_y = 0.0, coor_z = 0.0;

    element_s -> get_R(qua, R);
    
    for(int ii=0; ii<snLocBas; ++ii)
    {
      coor_x += sctrlPts_x[ii] * R[ii];
      coor_y += sctrlPts_y[ii] * R[ii];
      coor_z += sctrlPts_z[ii] * R[ii];

      sol_u  += solu[ii] * R[ii];
      sol_v  += solv[ii] * R[ii];
      sol_w  += solw[ii] * R[ii];
    }

    exact_wss( coor_x, coor_y, coor_z, t, exa_u, exa_v, exa_w );
    
    error += (sol_u - exa_u) * (sol_u - exa_u) * gwts;
    error += (sol_v - exa_v) * (sol_v - exa_v) * gwts;
    error += (sol_w - exa_w) * (sol_w - exa_w) * gwts;
  }

  return error;
}
