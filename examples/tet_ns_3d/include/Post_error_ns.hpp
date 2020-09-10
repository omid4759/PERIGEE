#ifndef POST_ERROR_NS_HPP
#define POST_ERROR_NS_HPP

#include <complex_bessel.h>
#include "Math_Tools.hpp"
#include "Matrix_3x3.hpp"
#include "FEAElement.hpp"

namespace POST_T_NS
{
  // Define some constants
  const double mu = 4.0e-2;
  const double rho0 = 1.0;

  const double R_pipe = 0.3;                                     // pipe radius
  const double omega  = MATH_T::PI * 2.0 / 1.1;                  // freqency
  const std::complex<double> i1(0.0, 1.0);
  const std::complex<double> i1_0d5( 0.707106781186547, 0.707106781186547);
  const std::complex<double> i1_1d5(-0.707106781186547, 0.707106781186547);
  const auto Omega    = std::sqrt(rho0 * omega / mu) * R_pipe;   // womersley number 
  const auto Lambda   = i1_1d5 * Omega;

  // Define pressure Fourier coefficients
  const double k0 = -21.0469;
  const std::complex<double> k1( -33.0102, 42.9332 );

  const auto coef1 = i1 * k1 / (rho0 * omega);

  // Exact solutions and their gradients
  double exact_pres( const double &x, const double &y, const double &z,
      const double &t );

  void exact_velo( const double &x, const double &y, const double &z,
      const double &t, double &val_x, double &val_y, double &val_z );

  void exact_grad_pres( const double &x, const double &y, const double &z,
      const double &t, double &val_x, double &val_y, double &val_z );

  void exact_grad_velo( const double &x, const double &y, const double &z,
      const double &t, Matrix_3x3 &grad_velo );

  double exact_wss( const double &x, const double &y, const double &z,
      const double &t );

  double get_pres_l2_error( const double * const &sol,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      double * const &R,
      const double &t );

  double get_pres_h1_error( const double * const &sol,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      double * const &R,
      double * const &Rx,
      double * const &Ry,
      double * const &Rz,
      const double &t );

  double get_velo_l2_error( const double * const &solu,
      const double * const &solv, const double * const &solw,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      double * const &R,
      const double &t );

  double get_velo_h1_error( const double * const &solu,
      const double * const &solv, const double * const &solw,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      double * const &R,
      double * const &Rx,
      double * const &Ry,
      double * const &Rz,
      const double &t );

  double get_wss_l2_error( const double * const &solu,
      const double * const &solv, const double * const &solw,
      const FEAElement * const &element,
      const double * const &ectrlPts_x,
      const double * const &ectrlPts_y,
      const double * const &ectrlPts_z,
      const IQuadPts * const &quad,
      double * const &R,
      const double &t );
}

#endif
