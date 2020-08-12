#ifndef IPLOCASSEM_2X2BLOCK_HPP
#define IPLOCASSEM_2X2BLOCK_HPP
// ==================================================================
// IPLocAssem_2x2Block.hpp
// Interface for parallel local assembly routine specialized for 
// 2-by-2, four block system.
//                 A = [ A00, A01;
//                       A10, A11 ]
// If the first variable has n0 degree-of-freedoms, the second variable
// has n1 degree-of-freedoms, then the size of the sub matrices should
// be: 
// A00 : n0 x n0
// A01 : n0 x n1
// A10 : n1 x n0
// A11 : n1 x n1.
//
// The Residual vector will be 
//                 R = [ R0;
//                       R1 ]
// R0 : n0,
// R1 : n1.
//
// All sub-matrices and sub-vectors are stored in the array vector
// of PetscScalar type.
// 
// If the input solution vectors are vec_0 and vec_1, they typically
// represent pressure and displacement/velocity.
// If the input solution vectors are vec_0, vec_1, vec_2, they
// typically represent displacement, pressure, and velocity.
//
// Author: Ju Liu
// Date: Jan 17 2018.
// ==================================================================
#include "FEAElement.hpp"
#include "IQuadPts.hpp"

class IPLocAssem_2x2Block
{
  public:
    IPLocAssem_2x2Block()
    {
      Tangent00 = nullptr;
      Tangent01 = nullptr;
      Tangent10 = nullptr;
      Tangent11 = nullptr;

      sur_Tangent00 = nullptr;
      sur_Tangent01 = nullptr;
      sur_Tangent10 = nullptr;
      sur_Tangent11 = nullptr;

      Residual0 = nullptr;
      Residual1 = nullptr;
      
      sur_Residual0 = nullptr;
      sur_Residual1 = nullptr;
    };

    virtual ~IPLocAssem_2x2Block(){};

    PetscScalar * Tangent00; // A00
    PetscScalar * Tangent01; // A01
    PetscScalar * Tangent10; // A10
    PetscScalar * Tangent11; // A11

    PetscScalar * sur_Tangent00; // sur_A00
    PetscScalar * sur_Tangent01; // sur_A01
    PetscScalar * sur_Tangent10; // sur_A10
    PetscScalar * sur_Tangent11; // sur_A11

    PetscScalar * Residual0; // R0
    PetscScalar * Residual1; // R1

    PetscScalar * sur_Residual0; // sur_R0
    PetscScalar * sur_Residual1; // sur_R1
    
    // -------------------------------------------------------------- 
    // Get the degree-of-freedom of this problem. For segregated 
    // algorithms this returns the fully-coupled multiphysics 
    // problem's dof.
    // -------------------------------------------------------------- 
    virtual int get_dof() const = 0;

    // -------------------------------------------------------------- 
    // Get the degree of freedom of the full matrix (i.e. n0 + n1).
    // In most problems, this value is the same as the get_dof value;
    // In the solid dynamics with kinematic segregated, this returns 4
    // (3 for velocity plus 1 for pressure); while get_dof returns 7 with
    // 3 additional displacement variables.
    // -------------------------------------------------------------- 
    virtual int get_dof_mat() const = 0;

    // -------------------------------------------------------------- 
    // Return the first variable's dof.
    // -------------------------------------------------------------- 
    virtual int get_dof_mat_0() const = 0;
    
    // -------------------------------------------------------------- 
    // Return the second variable's dof.
    // -------------------------------------------------------------- 
    virtual int get_dof_mat_1() const = 0;

    // -------------------------------------------------------------- 
    // Return the number of ebc functions implemented inside this local
    // assembly routine.
    // -------------------------------------------------------------- 
    virtual int get_num_ebc_fun() const
    {
      SYS_T::commPrint("Warning: IPLocAssem_2x2Block::get_num_ebc_fun is not implemented. \n");
      return 0;
    }
    
    // Assign all values in Tangent and Residual as 0.0.
    virtual void Zero_Tangent_Residual() = 0;

    virtual void Zero_sur_Tangent_Residual()
    {
      SYS_T::print_fatal("Error: Zero_sur_Tangent_Residual is not implemented.\n");
    }

    // Assign all values in Residual as 0.0.
    virtual void Zero_Residual() = 0;

    virtual void Zero_sur_Residual()
    {
      SYS_T::print_fatal("Error: Zero_sur_Residual is not implemented. \n");
    }

    // Generate nonzero pattern of all the sparse matrices.
    virtual void Assem_Estimate() = 0;

    // Assembly the Residuals: Residual0 and Residual1. 
    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &dot_vec,
        const double * const &vec,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: IPLocAssem_2x2Block::Assem_Residual is not implemented. \n");}

    virtual void Assem_Residual(
        const double &time, const double &dt,
        const double * const &dot_vec_0,
        const double * const &dot_vec_1,
        const double * const &dot_vec_2,
        const double * const &vec_0,
        const double * const &vec_1,
        const double * const &vec_2,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: IPLocAssem_2x2Block::Assem_Residual is not implemented. \n");}

    // Assembly the two Residuals and the four matrice blocks.
    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_vec,
        const double * const &vec,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: IPLocAssem_2x2Block::Assem_Tangent_Residual is not implemented. \n");}

    virtual void Assem_Tangent_Residual(
        const double &time, const double &dt,
        const double * const &dot_vec_0,
        const double * const &dot_vec_1,
        const double * const &dot_vec_2,
        const double * const &vec_0,
        const double * const &vec_1,
        const double * const &vec_2,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: IPLocAssem_2x2Block::Assem_Tangent_Residual is not implemented. \n");}

    // Assembly the two Residuals and the four matrice blocks for mass matrices.
    virtual void Assem_Mass_Residual(
        const double * const &vec,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: IPLocAssem_2x2Block::Assem_Mass_Residual is not implemented. \n");}

    virtual void Assem_Mass_Residual(
        const double * const &vec_0,
        const double * const &vec_1,
        const double * const &vec_2,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: IPLocAssem_2x2Block::Assem_Mass_Residual is not implemented. \n");}

    // Perform surface integration for elemental BC id ebc_id.
    virtual void Assem_Residual_EBC(
        const int &ebc_id,
        const double &time, const double &dt,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: IPLocAssem_2x2Block::Assem_Residual_EBC is not implemented.\n");}


    virtual double get_flowrate( const double * const &vec,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {
      SYS_T::commPrint("Warning: get_flowrate() is not implemented. \n");
      return 0.0;
    }

    virtual void get_pressure_area( const double * const &vec,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad,
        double &pres, double &area )
    {
      SYS_T::commPrint("Warning: get_pressure_area() is not implemented. \n");
    }

    virtual void Assem_Residual_EBC_Resistance(
        const int &ebc_id,
        const double &val,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Residual_EBC_Resistance is not implemented.\n");}

    virtual void Assem_Residual_BackFlowStab(
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Residual_BackFlowStab is not implemented.\n");}

    virtual void Assem_Tangent_Residual_BackFlowStab(
        const double &dt,
        const double * const &dot_sol,
        const double * const &sol,
        FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const IQuadPts * const &quad )
    {SYS_T::commPrint("Warning: this Assem_Tangent_Residual_BackFlowStab is not implemented.\n");}


    virtual double get_model_para_1() const
    {
      SYS_T::commPrint("Warning: get_model_para_1() is not implemented. \n");
      return 0.0;
    }

    virtual double get_model_para_2() const
    {
      SYS_T::commPrint("Warning: get_model_para_2() is not implemented. \n");
      return 0.0;
    }

};

#endif
