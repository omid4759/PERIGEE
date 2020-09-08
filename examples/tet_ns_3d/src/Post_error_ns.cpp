#include "Post_error_ns.hpp"

double POST_T_NS::exact_pres( const double &x, const double &y, const double &z,
          const double &t )
{
  return k0 * z + std::real( k1 * z * exp(i1*omega*t) ); 
}


void POST_T_NS::exact_velo( const double &x, const double &y, const double &z,
    const double &t, double &val_x, double &val_y, double &val_z )
{
  const double r     = sqrt(x*x+y*y);    // radial coord
  const auto   xi    = Lambda * r / R;
  const auto bes_top = sp_bessel::besselJ(0, xi);
  const auto bes_bot = sp_bessel::besselJ(0, Lambda);

  val_x = 0.0;
  val_y = 0.0;
  val_z = k0*(x*x + y*y - R*R) / (4.0*mu) + std::real( coef1 * exp(i1*omega*t) * (1.0 - bes_top / bes_bot ) );;
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
  const double r     = sqrt(x*x+y*y);    // radial coord
  const auto   xi    = Lambda * r / R;
  const auto bes_top = sp_bessel::besselJ(1, xi);
  const auto bes_bot = sp_bessel::besselJ(0, Lambda);

  grad_velo(0,0) = 0.0;
  grad_velo(0,1) = 0.0;
  grad_velo(0,2) = 0.0;

  grad_velo(1,0) = 0.0;
  grad_velo(1,1) = 0.0;
  grad_velo(1,2) = 0.0;

  grad_velo(2,0) = k0 * x / (2.0*mu) + std::real( coef1 * exp(i1*omega*t) * bes_top * i1_1d5 * Omega * x / (bes_bot * r * R) );
  grad_velo(2,1) = k0 * y / (2.0*mu) + std::real( coef1 * exp(i1*omega*t) * bes_top * i1_1d5 * Omega * y / (bes_bot * r * R) );
  grad_velo(2,2) = 0.0;
}

