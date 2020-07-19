#include "PGAssem_NLHeat_GenAlpha.hpp"

PGAssem_NLHeat_GenAlpha::PGAssem_NLHeat_GenAlpha( const IPLocAssem * const &locassem_ptr,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const APart_Node * const &pnode_ptr,
    const int &petsc_version_type )
{
  dof = locassem_ptr->get_dof();
  nLocBas = agmi_ptr->get_nLocBas();

  row_index = new PetscInt [dof * nLocBas];
  col_index = new PetscInt [dof * nLocBas];

  int sdegree = agmi_ptr->get_xdegree();
  int tdegree = agmi_ptr->get_ydegree();
  int udegree = agmi_ptr->get_zdegree();
  int nlocalnode = pnode_ptr->get_nlocalnode();

  int nz_prow = dof * (2*sdegree+1) * (2*tdegree+1) * (2*udegree+1);
  int nlocrow = dof * nlocalnode;

  switch(petsc_version_type)
  {
    case 0:
      Init_petsc_32(nz_prow, nlocrow);
      SYS_T::commPrint("===> PETSc-3.5.3: MatCreateAIJ called. \n");
      break;
    default:
      SYS_T::commPrint("Error: given petsc_version_type not implemented. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
      break;
  }

  // allocate the frequently used arrays in global assembly
  int nlgn = pnode_ptr->get_nlocghonode(); //number of local ghost nodes
  array_a = new double [nlgn * dof];
  array_b = new double [nlgn * dof];
  array_c = new double [nlgn * dof];
  array_d = new double [nlgn * dof];  

  local_a = new double [dof * nLocBas];
  local_b = new double [dof * nLocBas];
  local_c = new double [dof * nLocBas];
  local_d = new double [dof * nLocBas];  

  IEN_e = new int [nLocBas];

  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];
}


PGAssem_NLHeat_GenAlpha::PGAssem_NLHeat_GenAlpha( const IPLocAssem * const &locassem_ptr,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const APart_Node * const &pnode_ptr )
{
  dof = locassem_ptr->get_dof();
  nLocBas = agmi_ptr->get_nLocBas();

  row_index = new PetscInt [dof * nLocBas];
  col_index = new PetscInt [dof * nLocBas];

  const int sdegree = agmi_ptr->get_xdegree();
  const int tdegree = agmi_ptr->get_ydegree();
  const int udegree = agmi_ptr->get_zdegree();
  const int nlocalnode = pnode_ptr->get_nlocalnode();

  const int dnz = int ( 1.2 * dof * (2*sdegree+1) * (2*tdegree+1) * (2*udegree+1) );
  const int onz = int ( 1.2 * dof * (2*sdegree+1) * (2*tdegree+1) * (2*udegree+1) );
  const int nlocrow = dof * nlocalnode;

  Init_petsc_35(dnz, onz, nlocrow);
  SYS_T::commPrint("===> PETSc-3.5.3: MatCreateAIJ called. \n");

  Release_nonzero_err_str();

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE. \n");
  
  // allocate the frequently used arrays in global assembly
  int nlgn = pnode_ptr->get_nlocghonode(); //number of local ghost nodes
  array_a = new double [nlgn * dof];
  array_b = new double [nlgn * dof];

  local_a = new double [dof * nLocBas];
  local_b = new double [dof * nLocBas];

  IEN_e = new int [nLocBas];

  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];
}


PGAssem_NLHeat_GenAlpha::PGAssem_NLHeat_GenAlpha( const IPLocAssem * const &locassem_ptr,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &aien_ptr,
    const APart_Node * const &pnode_ptr,
    const IALocal_BC * const &part_bc )
{
  dof = locassem_ptr->get_dof();
  nLocBas = agmi_ptr->get_nLocBas();

  const int nlocalnode = pnode_ptr->get_nlocalnode();
  const int nlocrow = dof * nlocalnode;
  const int nElem = alelem_ptr->get_nlocalele();

  int * dnnz = new int [nlocrow];
  int * onnz = new int [nlocrow];

  SYS_T::commPrint("===> Estimate sparse nonzero structure. \n");
  Get_dnz_onz(nElem, aien_ptr, pnode_ptr, part_bc, dnnz, onnz);
  
  Init_petsc_35(dnnz, onnz, nlocrow);

  delete [] dnnz; delete [] onnz;

  // PETSc 3.5.3 and 3.6.0 use the same function call for creating Mat object  
  SYS_T::commPrint("===> PETSc-3.6.0: MatCreateAIJ called. \n");

  Release_nonzero_err_str();

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE. \n");

  // allocate the frequently used arrays in global assembly
  int nlgn = pnode_ptr->get_nlocghonode();
  array_a = new double [nlgn * dof];
  array_b = new double [nlgn * dof];

  row_index = new PetscInt [dof * nLocBas];
  col_index = new PetscInt [dof * nLocBas];
  
  local_a = new double [dof * nLocBas];
  local_b = new double [dof * nLocBas];

  IEN_e = new int [nLocBas];

  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];
}


PGAssem_NLHeat_GenAlpha::~PGAssem_NLHeat_GenAlpha()
{
  VecDestroy(&G);
  MatDestroy(&K);
  delete [] row_index;
  delete [] col_index;
  delete [] array_a;
  delete [] array_b;
  delete [] array_c;
  delete [] array_d;  
  delete [] local_a;
  delete [] local_b;
  delete [] local_c;
  delete [] local_d;  
  delete [] IEN_e;
  delete [] ectrl_x;
  delete [] ectrl_y;
  delete [] ectrl_z;
}

void PGAssem_NLHeat_GenAlpha::Init_petsc_32(const int &nonzero_per_row, const int &num_loc_row)
{
  MatCreateAIJ(PETSC_COMM_WORLD, num_loc_row, num_loc_row, PETSC_DECIDE,
      PETSC_DECIDE, nonzero_per_row, PETSC_NULL, nonzero_per_row, PETSC_NULL, &K);

  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, num_loc_row, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}


void PGAssem_NLHeat_GenAlpha::Init_petsc_35(const int &nonzero_per_row, const int &num_loc_row)
{
  MatCreateAIJ(PETSC_COMM_WORLD, num_loc_row, num_loc_row, PETSC_DECIDE,
      PETSC_DECIDE, nonzero_per_row, PETSC_NULL, nonzero_per_row, PETSC_NULL, &K);

  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, num_loc_row, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}


void PGAssem_NLHeat_GenAlpha::Init_petsc_35(const int &dnz, const int &onz, const int &num_loc_row)
{
  MatCreateAIJ(PETSC_COMM_WORLD, num_loc_row, num_loc_row, PETSC_DECIDE,
      PETSC_DECIDE, dnz, PETSC_NULL, onz, PETSC_NULL, &K);

  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, num_loc_row, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}


void PGAssem_NLHeat_GenAlpha::Init_petsc_35(const PetscInt * const &dnz,
    const PetscInt * const &onz, const int &num_loc_row)
{
  MatCreateAIJ(PETSC_COMM_WORLD, num_loc_row, num_loc_row, PETSC_DECIDE,
      PETSC_DECIDE, 0, dnz, 0, onz, &K);

  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, num_loc_row, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}


void PGAssem_NLHeat_GenAlpha::EssBC_KG(const IALocal_BC * const &bc_part, const int &field)
{
  // Check if dirichlet nodes exists within this partition
  // NOTE: we use ADD_VALUES here for matrix assembly, where we assumes that
  // the matrix is assemblyed with LID, which does nothing to the essential
  // boundary nodes, i.e., the boundary nodes' rows are zero rows.
  int local_dir = bc_part->get_Num_LD(field);
  if(local_dir > 0)
  {
    int row, col;
    for(int i=0; i<local_dir; ++i)
    {
      row = bc_part->get_LDN(field, i) * dof + field;
      col = row;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
    }
  }

  // Check if periodic slave nodes exists in this partition
  int local_sla = bc_part->get_Num_LP(field);
  if(local_sla > 0)
  {
    int row, col;
    for(int i=0; i<local_sla; ++i)
    {
      row = bc_part->get_LPSN(field, i) * dof + field;
      col = bc_part->get_LPMN(field, i) * dof + field;
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

void PGAssem_NLHeat_GenAlpha::EssBC_G(const IALocal_BC * const &bc_part, const int &field)
{
  int local_dir = bc_part->get_Num_LD(field);
  int local_sla = bc_part->get_Num_LP(field);
  int row;
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      row = bc_part->get_LDN(field, ii) * dof + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      row = bc_part->get_LPSN(field, ii) * dof + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}


// // ! GetLocal is moved to hpp to make it inline.
// void PGAssem_NLHeat_GenAlpha::GetLocal(const double * const &array, const int * const &IEN,
//     double * const &local_array) const
// {
//   for(int ii=0; ii<nLocBas; ++ii)
//   {
//     for(int jj=0; jj<dof; ++jj)
//       local_array[ii*dof+jj] = array[IEN[ii]*dof + jj];
//   }
// }


void PGAssem_NLHeat_GenAlpha::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr, 
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const IALocal_BC * const &bc_part )
{
  int nElem = alelem_ptr->get_nlocalele();
  int loc_dof = dof * nLocBas;
  int loc_index, lrow_index; // lcol_index;

  // loop over elements and insert 1.0 to every possible slots
  lassem_ptr->Assem_Estimate();

  for(int e=0; e<nElem; ++e)
  {
    for(int i=0; i<nLocBas; ++i)
    {
      loc_index  = lien_ptr->get_LIEN(e, i);
      
      for(int m=0; m<dof; ++m)
      {
        lrow_index = bc_part->get_LID( m, loc_index );

        row_index[dof * i + m] = dof * lrow_index + m;
        col_index[dof * i + m] = dof * lrow_index + m;
      }
    }
    MatSetValues(K, loc_dof, row_index, loc_dof, col_index,
        lassem_ptr->Tangent, ADD_VALUES);
  }

  for( int fie=0; fie<dof; ++fie)
    EssBC_KG( bc_part, fie );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}




void PGAssem_NLHeat_GenAlpha::Get_dnz_onz( const int &nElem,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const IALocal_BC * const &bc_part,
    PetscInt * const &dnz, PetscInt * const &onz ) const
{
  const int nlocalnode = node_ptr->get_nlocalnode();
  const int nnode = node_ptr->get_nlocghonode();

  std::vector<int> numLocNode;
  VEC_T::read_int_h5("NumLocalNode", "/", "nln", numLocNode);

  std::vector<unsigned int> nlist;
  nlist.clear();
  nlist.resize(numLocNode.size()+1);
  nlist[0] = 0;
  for(unsigned int ii=1; ii<=numLocNode.size(); ++ii)
    nlist[ii] = nlist[ii-1] + numLocNode[ii-1];

  // numlocNode is not needed anymore. nlist will provide the nodal 
  // partition info.
  VEC_T::clean(numLocNode);

  //for(int ii=0; ii<dof*nlocalnode; ++ii)
  //{
  //  dnz[ii] = 0; // Initialization: the num of nz are zero in each row.
  //  onz[ii] = 0;
  //}

  // This vector stores each row's diagonal col index and off-diagonal col index
  std::vector<int> interfun_d, interfun_o;

  // This MPI vector stores globally collected diagonal and off-diagonal col
  // number
  Vec vdnz, vonz;

  VecCreateMPI(PETSC_COMM_WORLD, dof * nlocalnode, PETSC_DETERMINE, &vdnz);
  VecCreateMPI(PETSC_COMM_WORLD, dof * nlocalnode, PETSC_DETERMINE, &vonz);

  VecSet(vdnz, 0.0);
  VecSet(vonz, 0.0);

  int row, col, ien_index, part_id_row, part_id_col;

  std::vector<int> elem4node;

  // loop for each row 
  for( int ii=0; ii<nnode; ++ii )
  {
    elem4node.clear();
    for(int ee=0; ee<nElem; ++ee)
    {
      if( lien_ptr->isNode_in_Elem(ee, ii) )
        elem4node.push_back(ee);
    }

    for( int mm=0; mm<dof; ++mm )
    {
      row = bc_part->get_LID( mm, ii );
      interfun_d.clear();
      interfun_o.clear();

      // only calculate for nondirichlet node
      if(row >= 0)
      {
        // based on nlist, find the row node's corresponding processor id.
        part_id_row = Get_part_id(nlist, row);
        for(unsigned int ei=0; ei<elem4node.size(); ++ei)
        {
          int ee = elem4node[ei];
          for(int kk=0; kk<nLocBas; ++kk)
          {
            ien_index = lien_ptr->get_LIEN(ee,kk);
            for( int nn=0; nn<dof; ++nn )
            {
              col = bc_part->get_LID( nn, ien_index );

              if( col>=0 )
              {
                // based on nlist, find the processor that the col node 
                // belong to. If they belong to the same processor, add 
                // to diagonal, otherwise, add to off-diagonal.
                part_id_col = Get_part_id(nlist, col);

                if( part_id_row == part_id_col )
                  interfun_d.push_back(dof*col + nn);
                else
                  interfun_o.push_back(dof*col + nn);
              }
            }
          }
        }
        // now the interfun_d and interfun_o have cached all the col that is
        // nonzero for this row, we sort them and remove repeated col number.
        VEC_T::sort_unique_resize(interfun_d);
        VEC_T::sort_unique_resize(interfun_o);

        // Finish calculating for each row, the real row index is
        // dof * row + mm
        VecSetValue(vdnz, dof*row+mm, double(interfun_d.size()), ADD_VALUES);
        VecSetValue(vonz, dof*row+mm, double(interfun_o.size()), ADD_VALUES);
      }
    }
  }

  VecAssemblyBegin(vdnz);
  VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz);
  VecAssemblyEnd(vonz);

  // We need to handle the Dirichlet and Periodic Slave nodes
  for(int mm=0; mm<dof; ++mm)
  {
    int local_dir = bc_part->get_Num_LD(mm);
    if(local_dir > 0)
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        int row = bc_part->get_LDN(mm, ii) * dof + mm;
        VecSetValue(vdnz, row, 1.0, INSERT_VALUES);
        VecSetValue(vonz, row, 0.0, INSERT_VALUES);
      }
    }

    int local_sla = bc_part->get_Num_LP(mm);

    if(local_sla > 0)
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        int row = bc_part->get_LPSN(mm,ii) * dof + mm;
        VecSetValue(vdnz, row, 2.0, INSERT_VALUES);
        VecSetValue(vonz, row, 2.0, INSERT_VALUES);
      }
    }
  }

  VecAssemblyBegin(vdnz);
  VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz);
  VecAssemblyEnd(vonz);

  // Get the global vector size
  PetscInt vec_size;
  VecGetSize(vdnz, &vec_size);

  // Now we get the globally collected dnz and onz number
  PetscScalar * array_d;
  PetscScalar * array_o;

  VecGetArray(vdnz, &array_d);
  for(int ii=0; ii<dof*nlocalnode; ++ii)
  {
    dnz[ii] = int(array_d[ii]);
    // the estimator above may overestimate for periodic master nodes.
    // if the number of nonzeros is greater than the dof * nlocalnode
    // reduce it to full diagonal rows. Otherwise PETSc will throw an
    // error message.
    if(dnz[ii] > dof*nlocalnode)
      dnz[ii] = dof * nlocalnode;                                   
  }
  VecRestoreArray(vdnz, &array_d);

  const int max_onz = vec_size - dof * nlocalnode;

  VecGetArray(vonz, &array_o);
  for(int ii=0; ii<dof*nlocalnode; ++ii)
  {
    onz[ii] = int(array_o[ii]);

    if(onz[ii] > max_onz)
      onz[ii] = max_onz;
  }
  VecRestoreArray(vonz, &array_o);

  VecDestroy(&vdnz);
  VecDestroy(&vonz);
}




//void PGAssem_NLHeat_GenAlpha::Assem_tangent_residual(
//    const PDNSolution * const &sol_a,
//    const PDNSolution * const &sol_b,
//    const double &curr_time,
//    const double &dt,
//    const ALocal_Elem * const &alelem_ptr,
//    IPLocAssem * const &lassem_ptr, 
//    const ALocal_IEN * const &lien_ptr,
//    const APart_Node * const &node_ptr,
//    const FEANode * const &fnode_ptr,
//    const AInt_Weight * const &wei_ptr,
//    const std::vector<FEAElement*> &eptr_array,
//    const IALocal_BC * const &bc_part )
//{
//  int nElem = alelem_ptr->get_nlocalele();
//  int loc_dof = dof * nLocBas;
//  int loc_index, lrow_index; // lcol_index;
//
//  sol_a->GetLocalArray( array_a, node_ptr );
//  sol_b->GetLocalArray( array_b, node_ptr );
//
//  for( int ee=0; ee<nElem; ++ee )
//  {
//    if( eptr_array[ee]->is_sizeNonzero() )
//    {
//      lien_ptr->get_LIEN_e(ee, IEN_e);
//      GetLocal(array_a, IEN_e, local_a); 
//      GetLocal(array_b, IEN_e, local_b);
//
//      fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);
//
//      lassem_ptr->Assem_Tangent_Residual(curr_time, dt, local_a, local_b,
//          eptr_array[ee], ectrl_x, ectrl_y, ectrl_z, wei_ptr);
//
//      for(int i=0; i<nLocBas; ++i)
//      {
//        loc_index = IEN_e[i];
//        //lcol_index = node_ptr->get_local_to_global( loc_index );
//
//        for(int m=0; m<dof; ++m)
//        {
//          lrow_index = bc_part->get_LID(m, loc_index);
//
//          row_index[dof * i + m] = dof * lrow_index + m;
//          col_index[dof * i + m] = dof * lrow_index + m;
//        }
//      }
//
//      MatSetValues(K, loc_dof, row_index, loc_dof, col_index,
//          lassem_ptr->Tangent, ADD_VALUES);
//
//      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
//    }
//  }
//
//  VecAssemblyBegin(G);
//  VecAssemblyEnd(G);
//
//  for(int fie = 0; fie<dof; ++fie)
//    EssBC_KG( bc_part, fie);
//
//  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
//  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
//  VecAssemblyBegin(G);
//  VecAssemblyEnd(G);
//}

//assemble when history variables exist
void PGAssem_NLHeat_GenAlpha::Assem_tangent_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const PDNSolution * const &sol_c, //ionic current
    const PDNSolution * const &sol_d, //dphi_ionic
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr, 
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const AInt_Weight * const &wei_ptr,
    const std::vector<FEAElement*> &eptr_array,
    const IALocal_BC * const &bc_part )
{
  std::cout << "PGAssem_NLHeat_GenAlpha" <<std::endl;
  int nElem = alelem_ptr->get_nlocalele();
  int loc_dof = dof * nLocBas;
  int loc_index, lrow_index; // lcol_index;

  sol_a->GetLocalArray( array_a, node_ptr );
  sol_b->GetLocalArray( array_b, node_ptr );
  sol_c->GetLocalArray( array_c, node_ptr );
  sol_d->GetLocalArray( array_d, node_ptr );  

  for( int ee=0; ee<nElem; ++ee )
  {
    if( eptr_array[ee]->is_sizeNonzero() )
    {
      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a); 
      GetLocal(array_b, IEN_e, local_b);
      GetLocal(array_c, IEN_e, local_c); 
      GetLocal(array_d, IEN_e, local_d);      

      fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      lassem_ptr->Assem_Tangent_Residual(curr_time, dt, local_a, local_b,
					 local_c, local_d, eptr_array[ee],
					 ectrl_x, ectrl_y, ectrl_z, wei_ptr);

      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];
        //lcol_index = node_ptr->get_local_to_global( loc_index );

        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);

          row_index[dof * i + m] = dof * lrow_index + m;
          col_index[dof * i + m] = dof * lrow_index + m;
        }
      }

      MatSetValues(K, loc_dof, row_index, loc_dof, col_index,
          lassem_ptr->Tangent, ADD_VALUES);

      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int fie = 0; fie<dof; ++fie)
    EssBC_KG( bc_part, fie);

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

//void PGAssem_NLHeat_GenAlpha::Assem_tangent_residual(
//    const PDNSolution * const &sol_a,
//    const PDNSolution * const &sol_b,
//    const double &curr_time,
//    const double &dt,
//    const ALocal_Elem * const &alelem_ptr,
//    IPLocAssem * const &lassem_ptr, 
//    const ALocal_IEN * const &lien_ptr,
//    const APart_Node * const &node_ptr,
//    const FEANode * const &fnode_ptr,
//    const AInt_Weight * const &wei_ptr,
//    const IALocal_meshSize * const &mSize,
//    const BernsteinBasis_Array * const &bs,
//    const BernsteinBasis_Array * const &bt,
//    const BernsteinBasis_Array * const &bu,
//    const IAExtractor * const &extractor,
//    const IALocal_BC * const &bc_part )
//{
//  const int nElem = alelem_ptr->get_nlocalele();
//  const int loc_dof = dof * nLocBas;
//  int loc_index, lrow_index; // lcol_index;
//
//  sol_a->GetLocalArray( array_a, node_ptr );
//  sol_b->GetLocalArray( array_b, node_ptr );
//
//  double ehx, ehy, ehz; // element mesh size in each direction
//
//  std::vector<double> ext_x, ext_y, ext_z; // extraction operator in each direction
//
//  double * ectrl_w = new double [nLocBas]; // element's weights
//
//  for( int ee=0; ee<nElem; ++ee )
//  {
//    //rmeshsize = mSize->get_meshsize(ee);
//
//    if( mSize->get_meshsize(ee) > 0.0 )
//    {
//      // Obtian the mesh size in element ee
//      ehx = mSize->get_hx(ee);
//      ehy = mSize->get_hy(ee);
//      ehz = mSize->get_hz(ee);
//
//      // Obtain the extractor in element ee
//      extractor->get_EXT_x(ee, ext_x);
//      extractor->get_EXT_y(ee, ext_y);
//      extractor->get_EXT_z(ee, ext_z);
//
//      lien_ptr->get_LIEN_e(ee, IEN_e);
//      GetLocal(array_a, IEN_e, local_a); 
//      GetLocal(array_b, IEN_e, local_b);
//
//      // Obtain the control points gemoetrical information
//      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);
//
//      // Call local assembly routine to generate element matrix and vector
//      lassem_ptr->Assem_Tangent_Residual(curr_time, dt, local_a, local_b,
//          ee, ehx, ehy, ehz, bs, bt, bu, ectrl_x, ectrl_y, ectrl_z, ectrl_w,
//          &ext_x[0], &ext_y[0], &ext_z[0], wei_ptr);
//
//
//      for(int i=0; i<nLocBas; ++i)
//      {
//        loc_index = IEN_e[i];
//
//        for(int m=0; m<dof; ++m)
//        {
//          lrow_index = bc_part->get_LID(m, loc_index);
//
//          row_index[dof * i + m] = dof * lrow_index + m;
//          col_index[dof * i + m] = dof * lrow_index + m;
//        }
//      }
//
//      MatSetValues(K, loc_dof, row_index, loc_dof, col_index,
//          lassem_ptr->Tangent, ADD_VALUES);
//
//      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
//    }
//  }
//
//  // Loop over top boundary face integral
//  // Due to the design of the boundary element info, one has to specify the
//  // location of boundary integral for the dof with index 0, even if the
//  // boundary condition is applied for other dofs. In the global
//  // assembly routine, I will only check the 0-index in the local bc object,
//  // if it indicates that the bc element number is greater than 0, boundary
//  // integral will be performed. This might be a bad design. However, if we
//  // adopt weak-imposition of dirichlet bc in the future, this does not matter
//  // at all. We will eventaully perform boundary integral for all dof's on all
//  // faces for all boundary conditions. The difference of difference bc will
//  // appear only in the local assembly routine.
//  for(int ei=0; ei<bc_part->get_NumLE_Top(0); ++ei)
//  {
//    // get the element index whose top surface needs to perform integral
//    int ee = bc_part->get_LTop_Elem(0, ei);
//    if( mSize->get_meshsize(ee) > 0.0 )
//    {
//      ehx = mSize->get_hx(ee);
//      ehy = mSize->get_hy(ee);
//      ehz = mSize->get_hz(ee);
//
//      extractor->get_EXT_x(ee, ext_x);
//      extractor->get_EXT_y(ee, ext_y);
//      extractor->get_EXT_z(ee, ext_z);
//
//      lien_ptr->get_LIEN_e(ee, IEN_e);
//      GetLocal(array_a, IEN_e, local_a); 
//      GetLocal(array_b, IEN_e, local_b);
//
//      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);
//
//      lassem_ptr->Assem_Residual_TopFace( curr_time, dt, local_a, local_b,
//          ee, ehx, ehy, ehz, ectrl_x, ectrl_y, ectrl_z, ectrl_w,
//          &ext_x[0], &ext_y[0], &ext_z[0] );
//
//      for(int i=0; i<nLocBas; ++i)
//      {
//        loc_index = IEN_e[i];
//        for(int m=0; m<dof; ++m)
//        {
//          lrow_index = bc_part->get_LID(m, loc_index);
//          row_index[dof * i + m] = dof * lrow_index + m;
//        } 
//      }
//      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
//    }
//  }
//  
//  for(int ei=0; ei<bc_part->get_NumLE_Bot(0); ++ei)
//  {
//    // get the element index whose top surface needs to perform integral
//    int ee = bc_part->get_LBottom_Elem(0, ei);
//    if( mSize->get_meshsize(ee) > 0.0 )
//    {
//      ehx = mSize->get_hx(ee);
//      ehy = mSize->get_hy(ee);
//      ehz = mSize->get_hz(ee);
//
//      extractor->get_EXT_x(ee, ext_x);
//      extractor->get_EXT_y(ee, ext_y);
//      extractor->get_EXT_z(ee, ext_z);
//
//      lien_ptr->get_LIEN_e(ee, IEN_e);
//      GetLocal(array_a, IEN_e, local_a); 
//      GetLocal(array_b, IEN_e, local_b);
//
//      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);
//
//      lassem_ptr->Assem_Residual_BotFace( curr_time, dt, local_a, local_b,
//          ee, ehx, ehy, ehz, ectrl_x, ectrl_y, ectrl_z, ectrl_w,
//          &ext_x[0], &ext_y[0], &ext_z[0] );
//
//      for(int i=0; i<nLocBas; ++i)
//      {
//        loc_index = IEN_e[i];
//        for(int m=0; m<dof; ++m)
//        {
//          lrow_index = bc_part->get_LID(m, loc_index);
//          row_index[dof * i + m] = dof * lrow_index + m;
//        } 
//      }
//      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
//    }
//  }
//
//
//  // The Rest boundary integrals need to be completed
//
//  VecAssemblyBegin(G);
//  VecAssemblyEnd(G);
//
//  // Imposing essential boundary conditions
//  for(int fie = 0; fie<dof; ++fie)
//    EssBC_KG( bc_part, fie);
//
//  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
//  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
//  VecAssemblyBegin(G);
//  VecAssemblyEnd(G);
//
//  delete [] ectrl_w;
//}


//void PGAssem_NLHeat_GenAlpha::Assem_residual(
//    const PDNSolution * const &sol_a,
//    const PDNSolution * const &sol_b,
//    const double &curr_time,
//    const double &dt,
//    const ALocal_Elem * const &alelem_ptr,
//    IPLocAssem * const &lassem_ptr, 
//    const ALocal_IEN * const &lien_ptr,
//    const APart_Node * const &node_ptr,
//    const FEANode * const &fnode_ptr,
//    const AInt_Weight * const &wei_ptr,
//    const std::vector<FEAElement*> &eptr_array,
//    const IALocal_BC * const &bc_part )
//{
//  int nElem = alelem_ptr->get_nlocalele();
//  int loc_dof = dof * nLocBas;
//  int loc_index, lrow_index;
//
//  sol_a->GetLocalArray( array_a, node_ptr );
//  sol_b->GetLocalArray( array_b, node_ptr );
//
//  for( int ee=0; ee<nElem; ++ee )
//  {
//    if( eptr_array[ee]->is_sizeNonzero() )
//    {
//      lien_ptr->get_LIEN_e(ee, IEN_e);
//      GetLocal(array_a, IEN_e, local_a); 
//      GetLocal(array_b, IEN_e, local_b);
//
//      fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);
//
//      lassem_ptr->Assem_Residual(curr_time, dt, local_a, local_b,
//          eptr_array[ee], ectrl_x, ectrl_y, ectrl_z, wei_ptr);
//
//      for(int i=0; i<nLocBas; ++i)
//      {
//        loc_index = IEN_e[i];
//
//        for(int m=0; m<dof; ++m)
//        {
//          lrow_index = bc_part->get_LID(m, loc_index);
//          row_index[dof * i + m] = dof * lrow_index + m;
//        }
//      }
//      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
//    }
//  }
//
//  VecAssemblyBegin(G);
//  VecAssemblyEnd(G);
//}

//assemble when history variables exist
void PGAssem_NLHeat_GenAlpha::Assem_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &sol_b,
    const PDNSolution * const &sol_c,//ionic current
    const PDNSolution * const &sol_d,//dphi_ionic       
    const double &curr_time,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr, 
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const AInt_Weight * const &wei_ptr,
    const std::vector<FEAElement*> &eptr_array,
    const IALocal_BC * const &bc_part )
{
  std::cout << "PGAssem_NLHeat_GenAlpha"<<std::endl;
  int nElem = alelem_ptr->get_nlocalele();
  int loc_dof = dof * nLocBas;
  int loc_index, lrow_index;

  sol_a->GetLocalArray( array_a, node_ptr );
  sol_b->GetLocalArray( array_b, node_ptr );
  sol_c->GetLocalArray( array_c, node_ptr );
  sol_d->GetLocalArray( array_d, node_ptr );

  for( int ee=0; ee<nElem; ++ee )
  {
    if( eptr_array[ee]->is_sizeNonzero() )
    {
      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a); 
      GetLocal(array_b, IEN_e, local_b);
      GetLocal(array_c, IEN_e, local_c); 
      GetLocal(array_d, IEN_e, local_d);
      
      fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      lassem_ptr->Assem_Residual(curr_time, dt, local_a, local_b,
				 local_c, local_d, eptr_array[ee],
				 ectrl_x, ectrl_y, ectrl_z, wei_ptr);

      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];

        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);
          row_index[dof * i + m] = dof * lrow_index + m;
        }
      }
      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


//void PGAssem_NLHeat_GenAlpha::Assem_residual(
//    const PDNSolution * const &sol_a,
//    const PDNSolution * const &sol_b,
//    const double &curr_time,
//    const double &dt,
//    const ALocal_Elem * const &alelem_ptr,
//    IPLocAssem * const &lassem_ptr, 
//    const ALocal_IEN * const &lien_ptr,
//    const APart_Node * const &node_ptr,
//    const FEANode * const &fnode_ptr,
//    const AInt_Weight * const &wei_ptr,
//    const IALocal_meshSize * const &mSize,
//    const BernsteinBasis_Array * const &bs,
//    const BernsteinBasis_Array * const &bt,
//    const BernsteinBasis_Array * const &bu,
//    const IAExtractor * const &extractor,
//    const IALocal_BC * const &bc_part )
//{
//  const int nElem = alelem_ptr->get_nlocalele();
//  const int loc_dof = dof * nLocBas;
//  int loc_index, lrow_index;
//
//  sol_a->GetLocalArray( array_a, node_ptr );
//  sol_b->GetLocalArray( array_b, node_ptr );
//
//  double ehx, ehy, ehz; // element mesh size in each direction
//
//  std::vector<double> ext_x, ext_y, ext_z; // extraction operator in each direction
//
//  double * ectrl_w = new double [nLocBas]; // element's weights
//
//  for( int ee=0; ee<nElem; ++ee )
//  {
//    if( mSize->get_meshsize(ee) > 0.0 )
//    {
//      // Obtain the mesh size in element ee
//      ehx = mSize->get_hx(ee);
//      ehy = mSize->get_hy(ee);
//      ehz = mSize->get_hz(ee);
//
//      // Obtain the extraction operator in element ee
//      extractor->get_EXT_x(ee, ext_x);
//      extractor->get_EXT_y(ee, ext_y);
//      extractor->get_EXT_z(ee, ext_z);
//
//      lien_ptr->get_LIEN_e(ee, IEN_e);
//      GetLocal(array_a, IEN_e, local_a); 
//      GetLocal(array_b, IEN_e, local_b);
//
//      // Obtain the geometrical infomation
//      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);
//
//      // Call local assembly routine to generate element matrix and vector
//      lassem_ptr->Assem_Residual(curr_time, dt, local_a, local_b,
//          ee, ehx, ehy, ehz, bs, bt, bu, ectrl_x, ectrl_y, ectrl_z, ectrl_w,
//          &ext_x[0], &ext_y[0], &ext_z[0], wei_ptr);
//
//      for(int i=0; i<nLocBas; ++i)
//      {
//        loc_index = IEN_e[i];
//
//        for(int m=0; m<dof; ++m)
//        {
//          lrow_index = bc_part->get_LID(m, loc_index);
//          row_index[dof * i + m] = dof * lrow_index + m;
//        }
//      }
//      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
//    }
//  }
//  
//  for(int ei=0; ei<bc_part->get_NumLE_Top(0); ++ei)
//  {
//    // get the element index whose top surface needs to perform integral
//    int ee = bc_part->get_LTop_Elem(0, ei);
//    if( mSize->get_meshsize(ee) > 0.0 )
//    {
//      ehx = mSize->get_hx(ee);
//      ehy = mSize->get_hy(ee);
//      ehz = mSize->get_hz(ee);
//
//      extractor->get_EXT_x(ee, ext_x);
//      extractor->get_EXT_y(ee, ext_y);
//      extractor->get_EXT_z(ee, ext_z);
//
//      lien_ptr->get_LIEN_e(ee, IEN_e);
//      GetLocal(array_a, IEN_e, local_a); 
//      GetLocal(array_b, IEN_e, local_b);
//
//      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);
//
//      lassem_ptr->Assem_Residual_TopFace( curr_time, dt, local_a, local_b,
//          ee, ehx, ehy, ehz, ectrl_x, ectrl_y, ectrl_z, ectrl_w,
//          &ext_x[0], &ext_y[0], &ext_z[0] );
//
//      for(int i=0; i<nLocBas; ++i)
//      {
//        loc_index = IEN_e[i];
//        for(int m=0; m<dof; ++m)
//        {
//          lrow_index = bc_part->get_LID(m, loc_index);
//          row_index[dof * i + m] = dof * lrow_index + m;
//        } 
//      }
//      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
//    }
//  }
//  
//  for(int ei=0; ei<bc_part->get_NumLE_Bot(0); ++ei)
//  {
//    // get the element index whose top surface needs to perform integral
//    int ee = bc_part->get_LBottom_Elem(0, ei);
//    if( mSize->get_meshsize(ee) > 0.0 )
//    {
//      ehx = mSize->get_hx(ee);
//      ehy = mSize->get_hy(ee);
//      ehz = mSize->get_hz(ee);
//
//      extractor->get_EXT_x(ee, ext_x);
//      extractor->get_EXT_y(ee, ext_y);
//      extractor->get_EXT_z(ee, ext_z);
//
//      lien_ptr->get_LIEN_e(ee, IEN_e);
//      GetLocal(array_a, IEN_e, local_a); 
//      GetLocal(array_b, IEN_e, local_b);
//
//      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);
//
//      lassem_ptr->Assem_Residual_BotFace( curr_time, dt, local_a, local_b,
//          ee, ehx, ehy, ehz, ectrl_x, ectrl_y, ectrl_z, ectrl_w,
//          &ext_x[0], &ext_y[0], &ext_z[0] );
//
//      for(int i=0; i<nLocBas; ++i)
//      {
//        loc_index = IEN_e[i];
//        for(int m=0; m<dof; ++m)
//        {
//          lrow_index = bc_part->get_LID(m, loc_index);
//          row_index[dof * i + m] = dof * lrow_index + m;
//        } 
//      }
//      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
//    }
//  }
//
//
//  // The Rest boundary integrals need to be completed
//
//  VecAssemblyBegin(G);
//  VecAssemblyEnd(G);
//
//  delete [] ectrl_w;
//}


void PGAssem_NLHeat_GenAlpha::Assem_mass_residual(
    const PDNSolution * const &sol_a,
    const PDNSolution * const &hist,
    const PDNTimeStep * const &time_info,
    const IonicModel * const &ionicmodel_ptr,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const AInt_Weight * const &wei_ptr,
    const std::vector<FEAElement*> &eptr_array,
    const IALocal_BC * const &bc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof * nLocBas;
  int loc_index, lrow_index;

  //FIRST: get ionic currents from initial potentials and hist.
  PDNSolution Iion(*hist); 
  PDNSolution dPhi_Iion(*hist);
  
  int node_num; 
  node_num =hist -> GetSize();
  double r_new_tmp, r_old, dt, new_soln, dPhi_Iion_tmp, Iion_tmp;
  dt=time_info->get_step();
  for (int count{ 0 }; count < node_num; ++count)
    {
      r_old = hist->GetValue(count);
      new_soln = sol_a->GetValue (count);

      ionicmodel_ptr->
	material_routine(r_old, dt,new_soln,
			 Iion_tmp, dPhi_Iion_tmp,
			 r_new_tmp);

      //hist_alpha.SetValue(count, r_new_alpha_tmp);
      //use negative below, to be consistent with krishnamoorthi
      //2013 quadrature paper and goktepe 2009 paper. 
      Iion.SetValue(count, -Iion_tmp);
      dPhi_Iion.SetValue(count, -dPhi_Iion_tmp);
    }
    
  sol_a->GetLocalArray( array_a, node_ptr );
  Iion.GetLocalArray( array_b, node_ptr );
  dPhi_Iion.GetLocalArray( array_c, node_ptr );
  
  for(int ee=0; ee<nElem; ++ee)
  {
    if( eptr_array[ee]->is_sizeNonzero() )
    {
      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a);
      GetLocal(array_b, IEN_e, local_b);//iion
      GetLocal(array_c, IEN_e, local_c);//d_iion

      fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      lassem_ptr->Assem_Mass_Residual(local_a, local_b, local_c,
				      eptr_array[ee], ectrl_x, ectrl_y,
				      ectrl_z, wei_ptr); 

      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];

        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);

          row_index[dof * i + m] = dof * lrow_index + m;
          col_index[dof * i + m] = dof * lrow_index + m;
        }
      }
      MatSetValues(K, loc_dof, row_index, loc_dof, col_index,
          lassem_ptr->Tangent, ADD_VALUES);

      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int fie = 0; fie<dof; ++fie)
    EssBC_KG( bc_part, fie);

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


//void PGAssem_NLHeat_GenAlpha::Assem_mass_residual(
//    const PDNSolution * const &sol_a,
//    const ALocal_Elem * const &alelem_ptr,
//    IPLocAssem * const &lassem_ptr,
//    const ALocal_IEN * const &lien_ptr,
//    const APart_Node * const &node_ptr,
//    const FEANode * const &fnode_ptr,
//    const AInt_Weight * const &wei_ptr,
//    const IALocal_meshSize * const &mSize,
//    const BernsteinBasis_Array * const &bs,
//    const BernsteinBasis_Array * const &bt,
//    const BernsteinBasis_Array * const &bu,
//    const IAExtractor * const &extractor,
//    const IALocal_BC * const &bc_part )
//{
//  const int nElem = alelem_ptr->get_nlocalele();
//  const int loc_dof = dof * nLocBas;
//  int loc_index, lrow_index;
//
//  double ehx, ehy, ehz; // element mesh size in each direction
//
//  std::vector<double> ext_x, ext_y, ext_z; // extraction operator in each direction
//
//  double * ectrl_w = new double [nLocBas]; // element's weights
//
//  sol_a->GetLocalArray( array_a, node_ptr );
//
//  for(int ee=0; ee<nElem; ++ee)
//  {
//    if( mSize->get_meshsize(ee) > 0.0 )
//    {
//      // obtain the mesh size in element ee
//      ehx = mSize->get_hx(ee);
//      ehy = mSize->get_hy(ee);
//      ehz = mSize->get_hz(ee);
//
//      // obtain the extraction operator in element ee
//      extractor->get_EXT_x(ee, ext_x);
//      extractor->get_EXT_y(ee, ext_y);
//      extractor->get_EXT_z(ee, ext_z);
//
//      lien_ptr->get_LIEN_e(ee, IEN_e);
//      GetLocal(array_a, IEN_e, local_a);
//
//      // obtain geometrical info
//      fnode_ptr->get_ctrlPts_xyzw(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z, ectrl_w);
//
//      lassem_ptr->Assem_Mass_Residual(local_a, ee, ehx, ehy, ehz, bs, bt, bu, 
//          ectrl_x, ectrl_y, ectrl_z, ectrl_w, &ext_x[0], &ext_y[0], &ext_z[0], wei_ptr); 
//
//      for(int i=0; i<nLocBas; ++i)
//      {
//        loc_index = IEN_e[i];
//
//        for(int m=0; m<dof; ++m)
//        {
//          lrow_index = bc_part->get_LID(m, loc_index);
//
//          row_index[dof * i + m] = dof * lrow_index + m;
//          col_index[dof * i + m] = dof * lrow_index + m;
//        }
//      }
//      MatSetValues(K, loc_dof, row_index, loc_dof, col_index,
//          lassem_ptr->Tangent, ADD_VALUES);
//
//      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
//    }
//  }
//  VecAssemblyBegin(G);
//  VecAssemblyEnd(G);
//
//  for(int fie = 0; fie<dof; ++fie)
//    EssBC_KG( bc_part, fie);
//
//  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
//  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
//  VecAssemblyBegin(G);
//  VecAssemblyEnd(G);
//
//  delete [] ectrl_w;
//}




//EOF
