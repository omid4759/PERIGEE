#ifndef IMATERIALMODEL_HPP
#define IMATERIALMODEL_HPP
// ==================================================================
// IMaterialModel.hpp
// 
// Interface for electroactive material models .
// Aliev-Panfilov model from Goktepe&Kuhl,2009
//
// Date: May 25 2020
// Author: Oguz Ziya Tikenogullari, Ju Liu
// Contact: o.z.tikenogullari@gmail.com, liujuy@gmail.com
// ==================================================================
#include "Sys_Tools.hpp"
#include "Tensor4_3D.hpp"
#include "Math_Tools.hpp"

class IonicModel
{
public:
  IonicModel();
  
  virtual ~IonicModel();
  
  virtual void print_info() const; //const =0 ;

  virtual double get_diso() const;

  virtual double get_dani() const;

  virtual double get_chi() const;

  virtual double get_C_m() const;
  
  virtual void material_routine(const double &r_old,
				const double &dt,
				const double &Phi,
				double &f_Phi,
				double &dP_fP,
				double &r_new) const;
private:
  const double ap_1, ap_2, ap_3, m1, m2, alpha, gamma,
    b, c, d_iso, d_ani, tol, chi, C_m;

//  virtual void get_PK(const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S) = 0;
//
//  virtual void get_PK_Stiffness(const Matrix_3x3 &F, Matrix_3x3 &P, 
//				Matrix_3x3 &S, Tensor4_3D &CC) = 0;
//
//  // Input: F : deformation gradient
//  // Output: P : 1st PK
//  //         S : 2nd PK
//  //         AA : F_iK F_jL C_KILJ
//  virtual void get_PK_FFStiffness( const Matrix_3x3 &F, Matrix_3x3 &P, 
//				   Matrix_3x3 &S, Tensor4_3D &AA )
//  {
//    get_PK_Stiffness(F, P, S, AA);
//    AA.MatMult_1(F);
//    AA.MatMult_3(F);
//  }
//
//  // Input: F : deformation gradient
//  // Output: sigma : Cauchy stress tensor
//  virtual void get_Cauchy_stress( const Matrix_3x3 &F, Matrix_3x3 &sigma )
//  {
//    Matrix_3x3 P, S;
//    get_PK(F, P, S);
//    Matrix_3x3 Ft(F); Ft.transpose();
//    sigma.MatMult(P, Ft);
//    sigma.scale( (1.0/F.det()) );
//  }
//
//  // Input: F : deformation gradient
//  // Output: sigma : Cauchy stress tensor
//  //         aa : J^{-1} F_iI F_jJ F_kK F_lL C_IJKL
//  virtual void get_Cauchy_stiffness( const Matrix_3x3 &F, Matrix_3x3 &sigma,
//				     Tensor4_3D &aa )
//  {
//    const double invJ = 1.0 / F.det();
//    Matrix_3x3 P, S;
//    get_PK_Stiffness(F, P, S, aa);
//    Matrix_3x3 Ft(F); Ft.transpose();
//    sigma.MatMult(P, Ft);
//    sigma.scale( invJ );
//
//    // Push forward the stiffness tensor
//    aa.MatMult_1(F); aa.MatMult_2(F); aa.MatMult_3(F); aa.MatMult_4(F);
//    aa.scale( invJ );
//  }
//
//  virtual double get_strain_energy( const Matrix_3x3 &F )
//  {
//    SYS_T::commPrint("Warning: IonicModel is not implemented. \n");
//    return 0.0;
//  }

};

#endif
