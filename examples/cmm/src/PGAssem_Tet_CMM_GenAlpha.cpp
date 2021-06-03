#include "PGAssem_Tet_CMM_GenAlpha.hpp"

PGAssem_Tet_CMM_GenAlpha::PGAssem_Tet_CMM_GenAlpha(
    IPLocAssem * const &locassem_ptr,
    FEAElement * const &elements,
    const IQuadPts * const &quads,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &aien_ptr,
    const APart_Node * const &pnode_ptr,
    const ALocal_NodalBC * const &part_nbc,
    const ALocal_Ring_NodalBC * const &part_ringnbc,
    const ALocal_EBC * const &part_ebc,
    const IGenBC * const &gbc,
    const int &in_nz_estimate )
: nLocBas( agmi_ptr->get_nLocBas() ),
  dof_sol( pnode_ptr->get_dof() ),
  dof_mat( locassem_ptr->get_dof_mat() ),
  num_ebc( part_ebc->get_num_ebc() ),
  nlgn( pnode_ptr->get_nlocghonode() ),
  snLocBas( 0 ) 
{
  // Make sure the data structure is compatible
  SYS_T::print_fatal_if(dof_sol != locassem_ptr->get_dof(),
      "PGAssem_Tet_CMM_GenAlpha::dof_sol != locassem_ptr->get_dof(). \n");

  SYS_T::print_fatal_if(dof_mat != part_nbc->get_dofMat(),
      "PGAssem_Tet_CMM_GenAlpha::dof_mat != part_nbc->get_dofMat(). \n");

  // Make sure that the surface element's number of local basis are 
  // the same. This is an assumption in this assembly routine.
  if(num_ebc>0) snLocBas = part_ebc -> get_cell_nLocBas(0);
  
  for(int ebc_id=0; ebc_id < num_ebc; ++ebc_id){
    SYS_T::print_fatal_if(snLocBas != part_ebc->get_cell_nLocBas(ebc_id),
        "Error: in PGAssem_Tet_CMM_GenAlpha, snLocBas has to be uniform. \n");
  }

  const int nlocrow = dof_mat * pnode_ptr->get_nlocalnode();

  // Allocate the sparse matrix K
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, dof_mat*in_nz_estimate, NULL, dof_mat*in_nz_estimate, NULL, &K);

  // Allocate the vector G
  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, nlocrow, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE.\n");
  Release_nonzero_err_str();

  Assem_nonzero_estimate( alelem_ptr, locassem_ptr, 
      elements, quads, aien_ptr, pnode_ptr, part_nbc, part_ringnbc, part_ebc, gbc );

  // Obtain the precise dnz and onz count
  std::vector<int> Kdnz, Konz;
  PETSc_T::Get_dnz_onz(K, Kdnz, Konz);

  MatDestroy(&K); // Destroy the K with rough preallocation
 
  // Create Mat with precise preallocation 
  MatCreateAIJ(PETSC_COMM_WORLD, nlocrow, nlocrow, PETSC_DETERMINE,
      PETSC_DETERMINE, 0, &Kdnz[0], 0, &Konz[0], &K);
}


PGAssem_Tet_CMM_GenAlpha::~PGAssem_Tet_CMM_GenAlpha()
{
  VecDestroy(&G);
  MatDestroy(&K);
}


void PGAssem_Tet_CMM_GenAlpha::EssBC_KG(
    const ALocal_NodalBC * const &nbc_part,
    const int &field )
{
  const int local_dir = nbc_part->get_Num_LD(field);

  if(local_dir > 0)
  {
    for(int i=0; i<local_dir; ++i)
    {
      const int row = nbc_part->get_LDN(field, i) * dof_mat + field;
      
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, row, 1.0, ADD_VALUES);
    }
  }

  const int local_sla = nbc_part->get_Num_LPS(field);
  if(local_sla > 0)
  {
    for(int i=0; i<local_sla; ++i)
    {
      const int row = nbc_part->get_LPSN(field, i) * dof_mat + field;
      const int col = nbc_part->get_LPMN(field, i) * dof_mat + field;
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}


void PGAssem_Tet_CMM_GenAlpha::RingBC_KG(
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const int &dof, const int &nrow, const int &ncol,
    const PetscInt * const &row_index,
    const PetscInt * const &col_index,
    PetscScalar * const &Ke,
    PetscScalar * const &Ge )
{
  const int ringbc_type = ringnbc_part -> get_ringbc_type();

  // Clamped rings
  if( ringbc_type == 0 ) {}

  // Skew boundary conditions for in-plane motion of ring nodes
  else if( ringbc_type == 1 )
  {
    // Note: element tangent Ke from NatBC_Resis_KG isn't a square matrix,
    //       ncol >= nrow. Pad to square ncol x ncol.
    PetscScalar * rotmat_e = new PetscScalar [ncol * ncol] {};

    // Set diagonal entries to 1.0
    for( int ii = 0; ii < ncol; ++ii ) rotmat_e[ ii*ncol + ii ] = 1.0;
    
    bool ring_rows = false, ring_cols = false; 
    int pos = -1; 

    for( int ii = dof-1; ii < ncol; ii += dof )
    {
      // Use velo-Z dof to determine ring nodes
      const int dnode = ( col_index[ii] - 3 ) / dof_mat;
      if( ringnbc_part->is_inLDN( dnode, pos ) )
      {
        Matrix_3x3 Q = ringnbc_part->get_rotation_matrix( pos );
        Q.transpose(); // Skew-to-global transformation matrix

        // Only rotate velocity dofs
        for( int jj = 0; jj < 3; ++jj )
        {
          for( int kk = 0; kk < 3; ++kk )
            rotmat_e[ (ii-2+jj) * ncol + (ii-2+kk) ] = Q(jj, kk);
        }
 
        ring_cols = true;
      }
    }

    for( int ii = dof-1; ii < nrow; ii += dof )
    {
      // Use velo-Z dof to determine ring nodes
      const int dnode = ( row_index[ii] - 3 ) / dof_mat;
      if( ringnbc_part->is_inLDN( dnode, pos ) )
      {
        ring_rows = true;
        break;
      } 
    }    

    if( ring_rows && ring_cols )
    {
      // Rotate K: R^T * K * R = R_{ki} * K_{kl} * R_{lj}
      PetscScalar * Ke_temp = new PetscScalar [ncol * ncol] {};

      for( int ii = 0; ii < nrow * ncol; ++ii )
      {
        Ke_temp[ii] = Ke[ii];
        Ke[ii] = 0.0;
      }

      for( int ii = 0; ii < nrow; ++ii )
      {
        for( int jj = 0; jj < ncol; ++jj )
        {
          for( int kk = 0; kk < ncol; ++kk )
          {
            for(int ll = 0; ll < ncol; ++ll ) 
              Ke[ii*ncol+jj] += rotmat_e[kk*ncol+ii] * Ke_temp[kk*ncol+ll] * rotmat_e[ll*ncol+jj];
          }
        }
      }

      // Rotate G: R^T * G = R_{ji} * G_{j}
      PetscScalar * Ge_temp = new PetscScalar [ncol] {};
      for( int ii = 0; ii < nrow; ++ii )
      {
        Ge_temp[ii] = Ge[ii];
        Ge[ii] = 0.0;
      }

      for(int ii = 0; ii < nrow; ++ii )
      {
        for(int jj = 0; jj < ncol; ++jj )
          Ge[ii] += rotmat_e[jj*ncol+ii] * Ge_temp[jj];
      }

      delete [] Ke_temp; delete [] Ge_temp;
      Ke_temp = nullptr; Ge_temp = nullptr;
    }
    else if( ring_cols )
    {
      // Rotate K columns: K * R = K_{ik} * R_{kj}
      PetscScalar * Ke_temp = new PetscScalar [ncol * ncol] {};
      
      for( int ii = 0; ii < nrow * ncol; ++ii )
      {
        Ke_temp[ii] = Ke[ii];
        Ke[ii] = 0.0;
      }

      for( int ii = 0; ii < nrow; ++ii )
      {
        for( int jj = 0; jj < ncol; ++jj )
        {
          for( int kk = 0; kk < ncol; ++kk )
            Ke[ii*ncol+jj] += Ke_temp[ii*ncol+kk] * rotmat_e[kk*ncol+jj];
        }
      }

      delete [] Ke_temp; Ke_temp = nullptr;
    }

    delete [] rotmat_e; rotmat_e = nullptr;
  }
  else
    SYS_T::print_fatal("Error: this ringbc_type is not supported in PGAssem_Tet_CMM_GenAlpha.\n");
}


void PGAssem_Tet_CMM_GenAlpha::RingBC_G(
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const int &dof, const int &nrow,
    const PetscInt * const &row_index,
    PetscScalar * const &Ge )
{
  const int ringbc_type = ringnbc_part -> get_ringbc_type();

  // Clamped rings
  if( ringbc_type == 0 ) {}

  // Skew boundary conditions for in-plane motion of ring nodes
  else if( ringbc_type == 1 )
  {
    PetscScalar * rotmat_e = new PetscScalar [nrow * nrow] {};

    // Set diagonal entries to 1.0
    for( int ii = 0; ii < nrow; ++ii ) rotmat_e[ ii*nrow + ii ] = 1.0;  

    bool ring_rows = false;
    int pos = -1; 

    for( int ii = dof-1; ii < nrow; ii += dof )
    {
      // Use velo-Z dof to determine ring nodes
      const int dnode = ( row_index[ii] - 3 ) / dof_mat;
      if( ringnbc_part->is_inLDN( dnode, pos ) )
      {
        Matrix_3x3 Q = ringnbc_part->get_rotation_matrix( pos );
        Q.transpose(); // Skew-to-global transformation matrix

        // Only rotate velocity dofs
        for( int jj = 0; jj < 3; ++jj )
        {
          for( int kk = 0; kk < 3; ++kk )
            rotmat_e[ (ii-2+jj) * nrow + (ii-2+kk) ] = Q(jj, kk);
        }

        ring_rows = true;
      }
    }

    if( ring_rows )
    {
      // Rotate G: R^T * G = R_{ji} * G_{j}
      PetscScalar * Ge_temp = new PetscScalar [nrow] {};
      for( int ii = 0; ii < nrow; ++ii )
      {
        Ge_temp[ii] = Ge[ii];
        Ge[ii] = 0.0;
      }

      for(int ii = 0; ii < nrow; ++ii )
      {
        for(int jj = 0; jj < nrow; ++jj )
          Ge[ii] += rotmat_e[jj*nrow+ii] * Ge_temp[jj];
      }

      delete [] Ge_temp; Ge_temp = nullptr;
    }

    delete [] rotmat_e; rotmat_e = nullptr;
  }
  else
    SYS_T::print_fatal("Error: this ringbc_type is not supported in PGAssem_Tet_CMM_GenAlpha.\n");
}


void PGAssem_Tet_CMM_GenAlpha::EssBC_G( const ALocal_NodalBC * const &nbc_part, 
    const int &field )
{
  const int local_dir = nbc_part->get_Num_LD(field);
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      const int row = nbc_part->get_LDN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  const int local_sla = nbc_part->get_Num_LPS(field);
  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      const int row = nbc_part->get_LPSN(field, ii) * dof_mat + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}


void PGAssem_Tet_CMM_GenAlpha::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elements,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;

  lassem_ptr->Assem_Estimate();

  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  for(int e=0; e<nElem; ++e)
  {
    for(int i=0; i<nLocBas; ++i)
    {
      const int loc_index  = lien_ptr->get_LIEN(e, i);

      for(int m=0; m<dof_mat; ++m)
        row_index[dof_mat * i + m] = dof_mat * nbc_part->get_LID( m, loc_index ) + m;
    }

    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        lassem_ptr->Tangent, ADD_VALUES);
  }

  delete [] row_index; row_index = nullptr;

  // ====== ISL TEST: REPLACE NATBC_RESIS_KG WITH NATBC_G ======
  // // Create a temporary zero solution vector to feed Natbc_Resis_KG
  // PDNSolution * temp = new PDNSolution_NS( node_ptr, 0, false );

  // // 0.1 is an (arbitrarily chosen) nonzero time step size feeding the NatBC_Resis_KG 
  // NatBC_Resis_KG(0.1, temp, temp, lassem_ptr, elements, quad_s, nbc_part, ringnbc_part, ebc_part, gbc );

  // delete temp;
  // ===========================================================

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii=0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_Tet_CMM_GenAlpha::Assem_mass_residual(
    const PDNSolution * const &sol_a,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;

  double * array_a = new double [nlgn * dof_sol];
  double * local_a = new double [nLocBas * dof_sol];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol_a->GetLocalArray( array_a );

  for(int ee=0; ee<nElem; ++ee)
  {
    lien_ptr->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Mass_Residual( local_a, elementv,
        ectrl_x, ectrl_y, ectrl_z, quad_v );

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat*ii+mm] = dof_mat * nbc_part -> get_LID(mm, IEN_e[ii]) + mm;
    }

    RingBC_KG( ringnbc_part, dof_mat, nLocBas * dof_mat, nLocBas * dof_mat,
        row_index, row_index, lassem_ptr->Tangent, lassem_ptr->Residual );

    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        lassem_ptr->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] array_a; array_a = nullptr;
  delete [] local_a; local_a = nullptr;
  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_Tet_CMM_GenAlpha::Assem_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const PDNSolution * const &sol_wall_disp,
    const PDNSolution * const &dot_sol_np1,
    const PDNSolution * const &sol_np1,
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    FEAElement * const &elementw,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_part,
    const ALocal_EBC * const &ebc_wall_part,
    const IGenBC * const &gbc )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;

  double * array_a = new double [nlgn * dof_sol];
  double * array_b = new double [nlgn * dof_sol];
  double * local_a = new double [nLocBas * dof_sol];
  double * local_b = new double [nLocBas * dof_sol];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol_a->GetLocalArray( array_a );
  sol_b->GetLocalArray( array_b );

  for( int ee=0; ee<nElem; ++ee )
  {
    lien_ptr->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Residual(curr_time, dt, local_a, local_b,
        elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat*ii+mm] = dof_mat * nbc_part -> get_LID(mm, IEN_e[ii]) + mm;
    }

    RingBC_G( ringnbc_part, dof_mat, nLocBas * dof_mat, row_index, lassem_ptr->Residual );

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] array_a; array_a = nullptr;
  delete [] array_b; array_b = nullptr;
  delete [] local_a; local_a = nullptr;
  delete [] local_b; local_b = nullptr;
  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;

  // Backflow stabilization residual contribution
  BackFlow_G( sol_a, sol_b, lassem_ptr, elements, quad_s, nbc_part, ringnbc_part, ebc_part );

  // Residual contribution from the thin-walled linear membrane in CMM
  WallMembrane_G( curr_time, dt, sol_a, sol_b, sol_wall_disp, lassem_ptr, elementw, quad_s, nbc_part, ringnbc_part, ebc_wall_part );

  // ====== ISL TEST: REPLACE NATBC_RESIS_KG WITH NATBC_G ======
  NatBC_G( curr_time, dt, lassem_ptr, elements, quad_s, nbc_part, ringnbc_part, ebc_part );

  // // Resistance type boundary condition
  // NatBC_Resis_G( dot_sol_np1, sol_np1, lassem_ptr, elements, quad_s, nbc_part, ebc_part, gbc );
  // ===========================================================

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_G( nbc_part, ii );

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_Tet_CMM_GenAlpha::Assem_tangent_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const PDNSolution * const &sol_wall_disp,
    const PDNSolution * const &dot_sol_np1,
    const PDNSolution * const &sol_np1,
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    FEAElement * const &elementw,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_part,
    const ALocal_EBC * const &ebc_wall_part,
    const IGenBC * const &gbc )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof_mat * nLocBas;

  double * array_a = new double [nlgn * dof_sol];
  double * array_b = new double [nlgn * dof_sol];
  double * local_a = new double [nLocBas * dof_sol];
  double * local_b = new double [nLocBas * dof_sol];
  int * IEN_e = new int [nLocBas];
  double * ectrl_x = new double [nLocBas];
  double * ectrl_y = new double [nLocBas];
  double * ectrl_z = new double [nLocBas];
  PetscInt * row_index = new PetscInt [nLocBas * dof_mat];

  sol_a->GetLocalArray( array_a );
  sol_b->GetLocalArray( array_b );

  for(int ee=0; ee<nElem; ++ee)
  {
    lien_ptr->get_LIEN(ee, IEN_e);
    GetLocal(array_a, IEN_e, local_a);
    GetLocal(array_b, IEN_e, local_b);

    fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    lassem_ptr->Assem_Tangent_Residual(curr_time, dt, local_a, local_b,
        elementv, ectrl_x, ectrl_y, ectrl_z, quad_v);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        row_index[dof_mat*ii + mm] = dof_mat*nbc_part->get_LID(mm, IEN_e[ii])+mm;
    }

    RingBC_KG( ringnbc_part, dof_mat, nLocBas * dof_mat, nLocBas * dof_mat,
        row_index, row_index, lassem_ptr->Tangent, lassem_ptr->Residual );

    MatSetValues(K, loc_dof, row_index, loc_dof, row_index,
        lassem_ptr->Tangent, ADD_VALUES);

    VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
  }

  delete [] array_a; array_a = nullptr;
  delete [] array_b; array_b = nullptr;
  delete [] local_a; local_a = nullptr;
  delete [] local_b; local_b = nullptr;
  delete [] IEN_e; IEN_e = nullptr;
  delete [] ectrl_x; ectrl_x = nullptr;
  delete [] ectrl_y; ectrl_y = nullptr;
  delete [] ectrl_z; ectrl_z = nullptr;
  delete [] row_index; row_index = nullptr;

  // Backflow stabilization residual & tangent contribution
  BackFlow_KG( dt, sol_a, sol_b, lassem_ptr, elements, quad_s, nbc_part, ringnbc_part, ebc_part );

  // Residual & tangent contributions from the thin-walled linear membrane in CMM
  WallMembrane_KG( curr_time, dt, sol_a, sol_b, sol_wall_disp, lassem_ptr, elementw, quad_s, nbc_part, ringnbc_part, ebc_wall_part );

  // ====== ISL TEST: REPLACE NATBC_RESIS_KG WITH NATBC_G ======
  NatBC_G( curr_time, dt, lassem_ptr, elements, quad_s, nbc_part, ringnbc_part, ebc_part );

  // // Resistance type boundary condition
  // NatBC_Resis_KG( dt, dot_sol_np1, sol_np1, lassem_ptr, elements, quad_s, nbc_part, ringnbc_part, ebc_part, gbc );
  // ===========================================================

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int ii = 0; ii<dof_mat; ++ii) EssBC_KG( nbc_part, ii );
  
  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_Tet_CMM_GenAlpha::NatBC_G( const double &curr_time, const double &dt,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_part )
{
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  PetscInt * srow_index = new PetscInt [dof_mat * snLocBas];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      lassem_ptr->Assem_Residual_EBC(ebc_id, curr_time, dt,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
          srow_index[dof_mat * ii + mm] = dof_mat * nbc_part -> get_LID(mm, LSIEN[ii]) + mm;
      }

      RingBC_G( ringnbc_part, dof_mat, snLocBas * dof_mat, srow_index, lassem_ptr->Residual );

      VecSetValues(G, dof_mat*snLocBas, srow_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }

  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
  delete [] srow_index; srow_index = nullptr;
}


void PGAssem_Tet_CMM_GenAlpha::BackFlow_G( 
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_part )
{
  double * array_a = new double [nlgn * dof_sol];
  double * array_b = new double [nlgn * dof_sol];
  double * local_as = new double [dof_sol * snLocBas];
  double * local_bs = new double [dof_sol * snLocBas];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  PetscInt * srow_index = new PetscInt [dof_mat * snLocBas];

  dot_sol->GetLocalArray( array_a );
  sol->GetLocalArray( array_b );

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      GetLocal(array_a, LSIEN, snLocBas, local_as);
      GetLocal(array_b, LSIEN, snLocBas, local_bs);

      lassem_ptr->Assem_Residual_BackFlowStab( local_as, local_bs,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
          srow_index[dof_mat * ii + mm] = dof_mat * nbc_part -> get_LID(mm, LSIEN[ii]) + mm;
      }

      RingBC_G( ringnbc_part, dof_mat, dof_mat * snLocBas, srow_index, lassem_ptr->sur_Residual );

      VecSetValues(G, dof_mat*snLocBas, srow_index, lassem_ptr->sur_Residual, ADD_VALUES);
    }
  }

  delete [] array_a; array_a = nullptr;
  delete [] array_b; array_b = nullptr;
  delete [] local_as; local_as = nullptr;
  delete [] local_bs; local_bs = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
  delete [] srow_index; srow_index = nullptr;
}


void PGAssem_Tet_CMM_GenAlpha::BackFlow_KG( const double &dt,
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_part )
{
  double * array_a = new double [nlgn * dof_sol];
  double * array_b = new double [nlgn * dof_sol];
  double * local_as = new double [dof_sol * snLocBas];
  double * local_bs = new double [dof_sol * snLocBas];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];
  PetscInt * srow_index = new PetscInt [dof_mat * snLocBas];

  dot_sol->GetLocalArray( array_a );
  sol->GetLocalArray( array_b );

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);

      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      GetLocal(array_a, LSIEN, snLocBas, local_as);
      GetLocal(array_b, LSIEN, snLocBas, local_bs);

      lassem_ptr->Assem_Tangent_Residual_BackFlowStab( dt, local_as, local_bs,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        for(int mm=0; mm<dof_mat; ++mm)
          srow_index[dof_mat * ii + mm] = dof_mat * nbc_part -> get_LID(mm, LSIEN[ii]) + mm;
      }

      RingBC_KG( ringnbc_part, dof_mat, dof_mat * snLocBas, dof_mat * snLocBas,
          srow_index, srow_index, lassem_ptr->sur_Tangent, lassem_ptr->sur_Residual );

      MatSetValues(K, dof_mat*snLocBas, srow_index, dof_mat*snLocBas, srow_index,
          lassem_ptr->sur_Tangent, ADD_VALUES);

      VecSetValues(G, dof_mat*snLocBas, srow_index, lassem_ptr->sur_Residual, ADD_VALUES);
    }
  }

  delete [] array_a; array_a = nullptr;
  delete [] array_b; array_b = nullptr;
  delete [] local_as; local_as = nullptr;
  delete [] local_bs; local_bs = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
  delete [] srow_index; srow_index = nullptr;
}


void PGAssem_Tet_CMM_GenAlpha::WallMembrane_G(
    const double &curr_time,
    const double &dt, 
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    const PDNSolution * const &sol_wall_disp,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_w,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_wall_part )
{
  const int dof_disp = 3; 

  double * array_a    = new double [nlgn * dof_mat ];
  double * array_b    = new double [nlgn * dof_mat];
  double * array_c    = new double [nlgn * dof_disp];
  double * local_as   = new double [snLocBas * dof_mat ];
  double * local_bs   = new double [snLocBas * dof_mat ];
  double * local_cs   = new double [snLocBas * dof_disp];
  int    * LSIEN      = new    int [snLocBas];
  double * sctrl_x    = new double [snLocBas];
  double * sctrl_y    = new double [snLocBas];
  double * sctrl_z    = new double [snLocBas];
  double * sthickness = new double [snLocBas];
  double * syoungsmod = new double [snLocBas];
  double * quaprestress = new double [ 6 * quad_s->get_num_quadPts() ];
  PetscInt * srow_index = new PetscInt [dof_mat * snLocBas];

  dot_sol->GetLocalArray( array_a );
  sol->GetLocalArray( array_b );

  sol_wall_disp->GetLocalArray( array_c );

  // wall has only one surface per the assumption in wall ebc
  const int ebc_id = 0;
  const int num_sele = ebc_wall_part -> get_num_local_cell(ebc_id);

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc_wall_part -> get_SIEN(ebc_id, ee, LSIEN);
    ebc_wall_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);
    ebc_wall_part -> get_thickness(ee, sthickness  );
    ebc_wall_part -> get_youngsmod(ee, syoungsmod  );
    ebc_wall_part -> get_prestress(ee, quaprestress);

    GetLocal(array_a, LSIEN, snLocBas, local_as);
    GetLocal(array_b, LSIEN, snLocBas, local_bs);

    GetLocal(array_c, LSIEN, snLocBas, dof_disp, local_cs);

    lassem_ptr->Assem_Residual_EBC_Wall( curr_time, dt, local_as, local_bs, local_cs,
        element_w, sctrl_x, sctrl_y, sctrl_z, sthickness, syoungsmod, quaprestress, quad_s);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        srow_index[dof_mat * ii + mm] = dof_mat * nbc_part -> get_LID(mm, LSIEN[ii]) + mm;
    }

    RingBC_G( ringnbc_part, dof_mat, dof_mat * snLocBas, srow_index, lassem_ptr->sur_Residual );

    VecSetValues(G, dof_mat*snLocBas, srow_index, lassem_ptr->sur_Residual, ADD_VALUES);
  }

  delete [] array_a;  array_a  = nullptr;
  delete [] array_b;  array_b  = nullptr;
  delete [] array_c;  array_c  = nullptr;
  delete [] local_as; local_as = nullptr;
  delete [] local_bs; local_bs = nullptr;
  delete [] local_cs; local_cs = nullptr;
  delete [] LSIEN;    LSIEN    = nullptr;
  delete [] sctrl_x;  sctrl_x  = nullptr;
  delete [] sctrl_y;  sctrl_y  = nullptr;
  delete [] sctrl_z;  sctrl_z  = nullptr;
  delete [] sthickness;   sthickness   = nullptr;
  delete [] syoungsmod;   syoungsmod   = nullptr;
  delete [] quaprestress; quaprestress = nullptr;
  delete [] srow_index;   srow_index   = nullptr;
}


void PGAssem_Tet_CMM_GenAlpha::WallMembrane_KG(
    const double &curr_time,
    const double &dt, 
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    const PDNSolution * const &sol_wall_disp,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_w,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_wall_part )
{
  const int dof_disp = 3; 

  double * array_a    = new double [nlgn * dof_mat ];
  double * array_b    = new double [nlgn * dof_mat ];
  double * array_c    = new double [nlgn * dof_disp];
  double * local_as   = new double [snLocBas * dof_mat ];
  double * local_bs   = new double [snLocBas * dof_mat ];
  double * local_cs   = new double [snLocBas * dof_disp];
  int    * LSIEN      = new    int [snLocBas];
  double * sctrl_x    = new double [snLocBas];
  double * sctrl_y    = new double [snLocBas];
  double * sctrl_z    = new double [snLocBas];
  double * sthickness = new double [snLocBas];
  double * syoungsmod = new double [snLocBas];
  double * quaprestress = new double [ 6 * quad_s->get_num_quadPts() ];
  PetscInt * srow_index = new PetscInt [dof_mat * snLocBas];

  dot_sol->GetLocalArray( array_a );
  sol->GetLocalArray( array_b );

  sol_wall_disp->GetLocalArray( array_c );

  // wall has only one surface per the assumption in wall ebc
  const int ebc_id = 0;
  const int num_sele = ebc_wall_part -> get_num_local_cell(ebc_id);

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc_wall_part -> get_SIEN(ebc_id, ee, LSIEN);
    ebc_wall_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);
    ebc_wall_part -> get_thickness(ee, sthickness  );
    ebc_wall_part -> get_youngsmod(ee, syoungsmod  );
    ebc_wall_part -> get_prestress(ee, quaprestress);

    GetLocal(array_a, LSIEN, snLocBas, local_as);
    GetLocal(array_b, LSIEN, snLocBas, local_bs);

    GetLocal(array_c, LSIEN, snLocBas, dof_disp, local_cs);

    lassem_ptr->Assem_Tangent_Residual_EBC_Wall( curr_time, dt, local_as, local_bs, local_cs,
        element_w, sctrl_x, sctrl_y, sctrl_z, sthickness, syoungsmod, quaprestress, quad_s);

    for(int ii=0; ii<snLocBas; ++ii)
    {
      for(int mm=0; mm<dof_mat; ++mm)
        srow_index[dof_mat * ii + mm] = dof_mat * nbc_part -> get_LID(mm, LSIEN[ii]) + mm;
    }

    RingBC_KG( ringnbc_part, dof_mat, dof_mat * snLocBas, dof_mat * snLocBas,
        srow_index, srow_index, lassem_ptr->sur_Tangent, lassem_ptr->sur_Residual );

    MatSetValues(K, dof_mat*snLocBas, srow_index, dof_mat*snLocBas, srow_index,
        lassem_ptr->sur_Tangent, ADD_VALUES);

    VecSetValues(G, dof_mat*snLocBas, srow_index, lassem_ptr->sur_Residual, ADD_VALUES);
  }

  delete [] array_a;  array_a  = nullptr;
  delete [] array_b;  array_b  = nullptr;
  delete [] array_c;  array_c  = nullptr;
  delete [] local_as; local_as = nullptr;
  delete [] local_bs; local_bs = nullptr;
  delete [] local_cs; local_cs = nullptr;
  delete [] LSIEN;    LSIEN    = nullptr;
  delete [] sctrl_x;  sctrl_x  = nullptr;
  delete [] sctrl_y;  sctrl_y  = nullptr;
  delete [] sctrl_z;  sctrl_z  = nullptr;
  delete [] sthickness; sthickness = nullptr;
  delete [] syoungsmod; syoungsmod = nullptr;
  delete [] srow_index; srow_index = nullptr;
}


void PGAssem_Tet_CMM_GenAlpha::Update_Wall_Prestress(
    const PDNSolution * const &sol_wall_disp,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_w,
    const IQuadPts * const &quad_s,
    ALocal_EBC * const &ebc_wall_part )
{
  const int dof_disp = 3;
  const int face_nqp = quad_s -> get_num_quadPts();

  double * array_b      = new double [nlgn * dof_disp];
  double * local_bs     = new double [snLocBas * dof_disp];
  int    * LSIEN        = new    int [snLocBas];

  double * syoungsmod   = new double [snLocBas];
  double * quaprestress = new double [6 * face_nqp];

  std::vector<Matrix_3x3> sigma; sigma.resize( face_nqp );

  sol_wall_disp->GetLocalArray( array_b );

  // wall has only one surface per the assumption in wall ebc
  const int ebc_id = 0;
  const int num_sele = ebc_wall_part -> get_num_local_cell(ebc_id);

  for(int ee=0; ee<num_sele; ++ee)
  {
    ebc_wall_part -> get_SIEN(ebc_id, ee, LSIEN);
    ebc_wall_part -> get_youngsmod(ee, syoungsmod  );
    ebc_wall_part -> get_prestress(ee, quaprestress);

    GetLocal(array_b, LSIEN, snLocBas, dof_disp, local_bs);

    lassem_ptr->get_Wall_CauchyStress( local_bs, element_w, syoungsmod, quad_s, sigma ); 

    // update prestress in Voigt notation (comps 11, 22, 33, 23, 13, 12)
    for(int qua=0; qua<face_nqp; ++qua)
    {
      quaprestress[6*qua]   += sigma[qua].xx();
      quaprestress[6*qua+1] += sigma[qua].yy();
      quaprestress[6*qua+2] += sigma[qua].zz();
      quaprestress[6*qua+3] += sigma[qua].yz();
      quaprestress[6*qua+4] += sigma[qua].xz();
      quaprestress[6*qua+5] += sigma[qua].xy();
    }

    ebc_wall_part -> set_prestress(ee, quaprestress);
  }

  delete [] array_b;      array_b      = nullptr;
  delete [] local_bs;     local_bs     = nullptr;
  delete [] LSIEN;        LSIEN        = nullptr;
  delete [] syoungsmod;   syoungsmod   = nullptr;
  delete [] quaprestress; quaprestress = nullptr;
}


double PGAssem_Tet_CMM_GenAlpha::Assem_surface_flowrate(
    const PDNSolution * const &vec,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_EBC * const &ebc_part,
    const int &ebc_id )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  vec -> GetLocalArray( array );

  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    ebc_part -> get_SIEN( ebc_id, ee, LSIEN);

    // Obtain the control points coordinates
    ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    esum += lassem_ptr -> get_flowrate( local, element_s, sctrl_x,
        sctrl_y, sctrl_z, quad_s );
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;

  double sum = 0.0;
  MPI_Allreduce(&esum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum;
}


double PGAssem_Tet_CMM_GenAlpha::Assem_surface_flowrate(
    const PDNSolution * const &vec,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_Inflow_NodalBC * const &infbc_part )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  vec -> GetLocalArray( array );

  const int num_sele = infbc_part -> get_num_local_cell();

  double esum = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    infbc_part -> get_SIEN( ee, LSIEN);

    // Obtain the control points coordinates
    infbc_part -> get_ctrlPts_xyz( ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    esum += lassem_ptr -> get_flowrate( local, element_s, sctrl_x,
        sctrl_y, sctrl_z, quad_s );
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;

  double sum = 0.0;
  MPI_Allreduce(&esum, &sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum;
}


double PGAssem_Tet_CMM_GenAlpha::Assem_surface_ave_pressure(
    const PDNSolution * const &vec,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_EBC * const &ebc_part,
    const int &ebc_id )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  vec -> GetLocalArray( array );

  const int num_sele = ebc_part -> get_num_local_cell(ebc_id);

  double val_pres = 0.0, val_area = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    ebc_part -> get_SIEN( ebc_id, ee, LSIEN);

    // Obtain the control points coordinates
    ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    double ele_pres, ele_area;

    lassem_ptr-> get_pressure_area( local, element_s, sctrl_x, sctrl_y,
        sctrl_z, quad_s, ele_pres, ele_area);

    val_pres += ele_pres;
    val_area += ele_area;
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;

  // Summation over CPUs
  double sum_pres = 0.0, sum_area = 0.0;

  MPI_Allreduce(&val_pres, &sum_pres, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&val_area, &sum_area, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum_pres / sum_area;
}


double PGAssem_Tet_CMM_GenAlpha::Assem_surface_ave_pressure(
    const PDNSolution * const &vec,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_Inflow_NodalBC * const &infbc_part )
{
  double * array = new double [nlgn * dof_sol];
  double * local = new double [snLocBas * dof_sol];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  vec -> GetLocalArray( array );

  const int num_sele = infbc_part -> get_num_local_cell();

  double val_pres = 0.0, val_area = 0.0;

  for(int ee=0; ee<num_sele; ++ee)
  {
    // Obtain the LSIEN array
    infbc_part -> get_SIEN( ee, LSIEN);

    // Obtain the control points coordinates
    infbc_part -> get_ctrlPts_xyz( ee, sctrl_x, sctrl_y, sctrl_z);

    // Obtain the solution vector in this element
    GetLocal(array, LSIEN, snLocBas, local);

    double ele_pres, ele_area;

    lassem_ptr-> get_pressure_area( local, element_s, sctrl_x, sctrl_y,
        sctrl_z, quad_s, ele_pres, ele_area);

    val_pres += ele_pres;
    val_area += ele_area;
  }

  delete [] array; array = nullptr;
  delete [] local; local = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;

  // Summation over CPUs
  double sum_pres = 0.0, sum_area = 0.0;

  MPI_Allreduce(&val_pres, &sum_pres, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&val_area, &sum_area, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return sum_pres / sum_area;
}


void PGAssem_Tet_CMM_GenAlpha::NatBC_Resis_G(
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  PetscScalar * Res = new PetscScalar [snLocBas * 3];
  PetscInt * srow_idx = new PetscInt [snLocBas * 3];
  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    // Calculate dot flow rate for face with ebc_id from solution vector dot_sol
    const double dot_flrate = Assem_surface_flowrate( dot_sol, lassem_ptr, 
        element_s, quad_s, ebc_part, ebc_id ); 

    // Calculate flow rate for face with ebc_id from solution vector sol
    const double flrate = Assem_surface_flowrate( sol, lassem_ptr,
        element_s, quad_s, ebc_part, ebc_id );

    // Get the (pressure) value on the outlet surface for traction evaluation    
    const double P_n   = gbc -> get_P0( ebc_id );
    const double P_np1 = gbc -> get_P( ebc_id, dot_flrate, flrate );

    // P_n+alpha_f
    // lassem_ptr->get_model_para_1() gives alpha_f 
    const double val = P_n + lassem_ptr->get_model_para_1() * (P_np1 - P_n);

    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);
    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      // Here, val is Pressure, and is used as the surface traction h = P I 
      // to calculate the boundary integral
      lassem_ptr->Assem_Residual_EBC_Resistance(ebc_id, val,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      for(int ii=0; ii<snLocBas; ++ii)
      {
        Res[3*ii+0] = lassem_ptr->Residual[4*ii+1];
        Res[3*ii+1] = lassem_ptr->Residual[4*ii+2];
        Res[3*ii+2] = lassem_ptr->Residual[4*ii+3];

        srow_idx[3*ii+0] = dof_mat * nbc_part->get_LID(1, LSIEN[ii]) + 1;
        srow_idx[3*ii+1] = dof_mat * nbc_part->get_LID(2, LSIEN[ii]) + 2;
        srow_idx[3*ii+2] = dof_mat * nbc_part->get_LID(3, LSIEN[ii]) + 3;
      }

      RingBC_G( ringnbc_part, 3, snLocBas * 3, srow_idx, Res );

      VecSetValues(G, snLocBas*3, srow_idx, Res, ADD_VALUES);
    }
  }

  delete [] Res; Res = nullptr; delete [] srow_idx; srow_idx = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
}


void PGAssem_Tet_CMM_GenAlpha::NatBC_Resis_KG(
    const double &dt,
    const PDNSolution * const &dot_sol,
    const PDNSolution * const &sol,
    IPLocAssem * const &lassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Ring_NodalBC * const &ringnbc_part,
    const ALocal_EBC * const &ebc_part,
    const IGenBC * const &gbc )
{
  const double a_f = lassem_ptr -> get_model_para_1();

  // dd_dv = dt x alpha_f x gamma
  const double dd_dv = dt * a_f * lassem_ptr->get_model_para_2();

  // Allocate the vector to hold the residual on each surface element
  PetscScalar * Res = new PetscScalar [snLocBas * 3];
  PetscInt * srow_idx = new PetscInt [snLocBas * 3];
  PetscScalar * Tan;
  PetscInt * scol_idx;
  Vector_3 out_n;
  std::vector<double> intNB;
  std::vector<int> map_Bj;

  int * LSIEN = new int [snLocBas];
  double * sctrl_x = new double [snLocBas];
  double * sctrl_y = new double [snLocBas];
  double * sctrl_z = new double [snLocBas];

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    // Calculate dot flow rate for face with ebc_id and MPI_Allreduce them
    // Here, dot_sol is the solution at time step n+1 (not n+alpha_f!)
    const double dot_flrate = Assem_surface_flowrate( dot_sol, lassem_ptr, 
        element_s, quad_s, ebc_part, ebc_id ); 

    // Calculate flow rate for face with ebc_id and MPI_Allreduce them
    // Here, sol is the solution at time step n+1 (not n+alpha_f!)
    const double flrate = Assem_surface_flowrate( sol, lassem_ptr,
        element_s, quad_s, ebc_part, ebc_id );

    // Get the (pressure) value on the outlet surface for traction evaluation    
    const double P_n   = gbc -> get_P0( ebc_id );
    const double P_np1 = gbc -> get_P( ebc_id, dot_flrate, flrate );

    // P_n+alpha_f 
    const double resis_val = P_n + a_f * (P_np1 - P_n);

    // Get the (potentially approximated) m := dP/dQ
    const double m_val = gbc -> get_m( ebc_id, dot_flrate, flrate );

    // Get the (potentially approximated) n := dP/d(dot_Q)
    const double n_val = gbc -> get_n( ebc_id, dot_flrate, flrate );

    // Define alpha_f * n + alpha_f * gamma * dt * m
    // coef a^t a enters as the consistent tangent for the resistance-type bc
    const double coef = a_f * n_val + dd_dv * m_val;

    const int num_face_nodes = ebc_part -> get_num_face_nodes(ebc_id);
    if(num_face_nodes > 0)
    {
      Tan = new PetscScalar [snLocBas * 3 * num_face_nodes * 3];
      scol_idx = new PetscInt [num_face_nodes * 3];
      out_n  = ebc_part -> get_outvec( ebc_id );
      intNB  = ebc_part -> get_intNA( ebc_id );
      map_Bj = ebc_part -> get_LID( ebc_id );
    }
    else
    {
      Tan = nullptr;
      scol_idx = nullptr;
    }

    const int num_sele = ebc_part -> get_num_local_cell(ebc_id);
    for(int ee=0; ee<num_sele; ++ee)
    {
      ebc_part -> get_SIEN(ebc_id, ee, LSIEN);
      ebc_part -> get_ctrlPts_xyz(ebc_id, ee, sctrl_x, sctrl_y, sctrl_z);

      // For here, we scale the int_NA nx/y/z by factor 1
      lassem_ptr->Assem_Residual_EBC_Resistance(ebc_id, 1.0,
          element_s, sctrl_x, sctrl_y, sctrl_z, quad_s);

      // Residual vector is scaled by the resistance value
      for(int ii=0; ii<snLocBas; ++ii)
      {
        Res[3*ii+0] = resis_val * lassem_ptr->Residual[4*ii+1];
        Res[3*ii+1] = resis_val * lassem_ptr->Residual[4*ii+2];
        Res[3*ii+2] = resis_val * lassem_ptr->Residual[4*ii+3];
      }

      for(int A=0; A<snLocBas; ++A)
      {
        for(int ii=0; ii<3; ++ii)
        {
          const int temp_row = (3*A+ii) * num_face_nodes * 3;
          for(int B=0; B<num_face_nodes; ++B)
          {
            // Residual[4*A+ii+1] is intNB[A]*out_n[ii]
            Tan[temp_row + 3*B + 0] = coef * lassem_ptr->Residual[4*A+ii+1] * intNB[B] * out_n.x();
            Tan[temp_row + 3*B + 1] = coef * lassem_ptr->Residual[4*A+ii+1] * intNB[B] * out_n.y();
            Tan[temp_row + 3*B + 2] = coef * lassem_ptr->Residual[4*A+ii+1] * intNB[B] * out_n.z();
          }
        }
      }

      for(int ii=0; ii<snLocBas; ++ii)
      {
        srow_idx[3*ii+0] = dof_mat * nbc_part->get_LID(1,LSIEN[ii]) + 1;
        srow_idx[3*ii+1] = dof_mat * nbc_part->get_LID(2,LSIEN[ii]) + 2;
        srow_idx[3*ii+2] = dof_mat * nbc_part->get_LID(3,LSIEN[ii]) + 3;
      }

      for(int ii=0; ii<num_face_nodes; ++ii)
      {
        scol_idx[ii*3+0] = dof_mat * map_Bj[ii*3+0] + 1;
        scol_idx[ii*3+1] = dof_mat * map_Bj[ii*3+1] + 2;
        scol_idx[ii*3+2] = dof_mat * map_Bj[ii*3+2] + 3;
      }

      RingBC_KG( ringnbc_part, 3, snLocBas * 3, num_face_nodes * 3, srow_idx, scol_idx, Tan, Res );

      MatSetValues(K, snLocBas*3, srow_idx, num_face_nodes*3, scol_idx, Tan, ADD_VALUES);
      VecSetValues(G, snLocBas*3, srow_idx, Res, ADD_VALUES);
    }

    if( num_face_nodes > 0 ) 
    {
      delete [] Tan; Tan = nullptr;
      delete [] scol_idx; scol_idx = nullptr;
    }
  }

  delete [] Res; Res = nullptr; delete [] srow_idx; srow_idx = nullptr;
  delete [] LSIEN; LSIEN = nullptr;
  delete [] sctrl_x; sctrl_x = nullptr;
  delete [] sctrl_y; sctrl_y = nullptr;
  delete [] sctrl_z; sctrl_z = nullptr;
}

void PGAssem_Tet_CMM_GenAlpha::print_2Darray( const double * const arr,
  const int &nrow, const int &ncol )
{
  for(int ii = 0; ii < nrow; ++ii)
  {
    for(int jj = 0; jj < ncol; ++jj)
    {
      std::cout << std::scientific << std::setprecision(3) << std::setw(10) << arr[ii * ncol + jj] << " ";
      // std::cout << arr[ii * ncol + jj] << " ";
    }
    std::cout << std::endl;
  }
}
// EOF
