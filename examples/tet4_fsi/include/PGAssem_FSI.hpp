#ifndef PGASSEM_FSI_HPP
#define PGASSEM_FSI_HPP
// ============================================================================
// PGAssem_FSI.hpp
// 
// Parallel global assembly for FSI problems using the unified continuum
// formulation and segregated algorithm.
// 
// Author: Ju Liu
// Date: Jan 2 2022
// ============================================================================
#include "IPGAssem.hpp"
#include "IPLocAssem_2x2Block.hpp"

class PGAssem_FSI : public IPGAssem
{
  public:
    PGAssem_FSI( 
        IPLocAssem_2x2Block * const &locassem_f_ptr,
        IPLocAssem_2x2Block * const &locassem_s_ptr,
        FEAElement * const &elements,
        const IQuadPts * const &quads,
        const IAGlobal_Mesh_Info * const &agmi_v,
        const IAGlobal_Mesh_Info * const &agmi_p,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_v,
        const ALocal_IEN * const &aien_p,
        const APart_Node * const &pnode_v,
        const APart_Node * const &pnode_p,
        const ALocal_NodalBC * const &part_nbc_v,
        const ALocal_NodalBC * const &part_nbc_p,
        const ALocal_EBC * const &part_ebc,
        const IGenBC * const &gbc,
        const int &in_nz_estimate = 60 );

    virtual ~PGAssem_FSI();


    // Assembly routine for the surface integrals for flow rates
    // and averaged pressure
    virtual double Assem_surface_flowrate(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_EBC * const &ebc_part,
        const int &ebc_id );
    
    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_EBC * const &ebc_part,
        const int &ebc_id );

    virtual double Assem_surface_flowrate(
        const PDNSolution * const &disp,
        const PDNSolution * const &velo,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_Inflow_NodalBC * const &infbc_part,
        const int &nbc_id );

    virtual double Assem_surface_ave_pressure(
        const PDNSolution * const &disp,
        const PDNSolution * const &pres,
        IPLocAssem_2x2Block * const &lassem_ptr,
        FEAElement * const &element_s,
        const IQuadPts * const &quad_s,
        const ALocal_Inflow_NodalBC * const &infbc_part,
        const int &nbc_id );
  
  private:
    const int nLocBas, snLocBas, num_ebc, nlgn_v, nlgn_p;

    void EssBC_KG( const ALocal_NodalBC * const &nbc_v, const ALocal_NodalBC * const &nbc_p );

    void EssBC_G( const ALocal_NodalBC * const &nbc_v, const ALocal_NodalBC * const &nbc_p );

    void GetLocal( const double * const &array, const int * const &IEN,
        const int &in_locbas, const int &in_dof, double * const &local_array) const
    {
      for(int ii=0; ii<in_locbas; ++ii)
        for(int jj=0; jj<in_dof; ++jj)
          local_array[ii * in_dof + jj] = array[IEN[ii] * in_dof + jj];
    }
};

#endif
