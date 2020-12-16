#include "PNonlinear_NS_Solver.hpp"

PNonlinear_NS_Solver::PNonlinear_NS_Solver(
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol,
    const int &input_max_iteration,
    const int &input_renew_freq,
    const int &input_renew_threshold )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq),
  nrenew_threshold(input_renew_threshold)
{
  // Generate the incremental solution vector used for update
  // the solution of the nonlinear algebraic system
  dot_step = new PDNSolution_NS( anode_ptr, 0, false );
}


PNonlinear_NS_Solver::~PNonlinear_NS_Solver()
{
  delete dot_step; dot_step = nullptr;
}


void PNonlinear_NS_Solver::print_info() const
{
  SYS_T::commPrint("----------------------------------------------------------- \n");
  SYS_T::commPrint("Nonlinear solver setted up:\n");
  SYS_T::commPrint("  relative tolerance: %e \n", nr_tol);
  SYS_T::commPrint("  absolute tolerance: %e \n", na_tol);
  SYS_T::commPrint("  divergence tolerance: %e \n", nd_tol);
  SYS_T::commPrint("  maximum iteration: %d \n", nmaxits);
  SYS_T::commPrint("  tangent matrix renew frequency: %d \n", nrenew_freq);
  SYS_T::commPrint("  tangent matrix renew threshold: %d \n", nrenew_threshold);
  SYS_T::commPrint("----------------------------------------------------------- \n");
}


void PNonlinear_NS_Solver::GenAlpha_Solve_NS(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &sol_base,
    const PDNSolution * const &pre_dot_sol,
    const PDNSolution * const &pre_sol,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ICVFlowRate * const flr_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &nbc_part,
    const ALocal_Inflow_NodalBC * const &infnbc_part,
    const ALocal_EBC * const &ebc_part,
    IGenBC * const &gbc,
    const Matrix_PETSc * const &bc_mat,
    FEAElement * const &elementv,
    FEAElement * const &elements,
    const IQuadPts * const &quad_v,
    const IQuadPts * const &quad_s,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PDNSolution * const &dot_sol,
    PDNSolution * const &sol,
    bool &conv_flag, int &nl_counter ) const
{
#ifdef PETSC_USE_LOG
  PetscLogEvent mat_assem_0_event, mat_assem_1_event;
  PetscLogEvent vec_assem_0_event, vec_assem_1_event;
  PetscLogEvent lin_solve_event;
  PetscClassId classid_assembly;
  PetscClassIdRegister("mat_vec_assembly", &classid_assembly);
  PetscLogEventRegister("assembly mat 0", classid_assembly, &mat_assem_0_event);
  PetscLogEventRegister("assembly mat 1", classid_assembly, &mat_assem_1_event);
  PetscLogEventRegister("assembly vec 0", classid_assembly, &vec_assem_0_event);
  PetscLogEventRegister("assembly vec 1", classid_assembly, &vec_assem_1_event);
  PetscLogEventRegister("lin_solve", classid_assembly, &lin_solve_event);
#endif

  // Initialize the counter and error
  nl_counter = 0;
  double residual_norm = 0.0, initial_norm = 0.0, relative_error = 0.0;

  // Gen-alpha parameters
  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

  // Same-Y predictor
  sol->Copy(*pre_sol);
  dot_sol->Copy(*pre_dot_sol);
  dot_sol->ScaleValue( (gamma-1.0)/gamma );

  // Define the dol_sol at alpha_m: dot_sol_alpha
  PDNSolution dot_sol_alpha(*pre_dot_sol);
  dot_sol_alpha.ScaleValue( 1.0 - alpha_m );
  dot_sol_alpha.PlusAX(*dot_sol, alpha_m);

  // Define the sol at alpha_f: sol_alpha
  PDNSolution sol_alpha(*pre_sol);
  sol_alpha.ScaleValue( 1.0 - alpha_f );
  sol_alpha.PlusAX( *sol, alpha_f );


  if(gbc-> UserLPM_Dirichlet_flag==false){
  // -------------------------------------------------
  // Update the inflow boundary values
    rescale_inflow_value(curr_time+dt, infnbc_part, flr_ptr, sol_base, sol);
    rescale_inflow_value(curr_time+alpha_f*dt, infnbc_part, flr_ptr, sol_base, &sol_alpha);

  }else{
  // call GenBC to update Dirichlet and Neumann boundary conditions
    get_gbc_pressure_flow(curr_time,dt,alpha_f,sol_base,sol,dot_sol,&sol_alpha,infnbc_part,ebc_part,gbc,lassem_ptr, gassem_ptr,elements,quad_s,new_tangent_flag);

  }
  // -------------------------------------------------

  // If new_tangent_flag == TRUE, update the tangent matrix;
  // otherwise, use the matrix from the previous time step
  if( new_tangent_flag )
  {
    gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(mat_assem_0_event, 0,0,0,0);
#endif

    gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part, gbc );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(mat_assem_0_event,0,0,0,0);
#endif

    SYS_T::commPrint("  --- M updated");

    // SetOperator will pass the tangent matrix to the linear solver and the
    // linear solver will generate the preconditioner based on the new matrix.
    lsolver_ptr->SetOperator( gassem_ptr->K );
  }
  else
  {
    gassem_ptr->Clear_G();

#ifdef PETSC_USE_LOG
    PetscLogEventBegin(vec_assem_0_event, 0,0,0,0);
#endif

    gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
        curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
        quad_v, quad_s, lien_ptr, anode_ptr,
        feanode_ptr, nbc_part, ebc_part, gbc );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(vec_assem_0_event,0,0,0,0);
#endif
  }

  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  SYS_T::commPrint("  Init res 2-norm: %e \n", initial_norm);

  // Now do consistent Newton-Raphson iteration
  do
  {
#ifdef PETSC_USE_LOG
    PetscLogEventBegin(lin_solve_event, 0,0,0,0);
#endif

    // solve the equation K dot_step = G
    lsolver_ptr->Solve( gassem_ptr->G, dot_step );

#ifdef PETSC_USE_LOG
    PetscLogEventEnd(lin_solve_event,0,0,0,0);
#endif

    bc_mat->MatMultSol( dot_step );

    nl_counter += 1;

    dot_sol->PlusAX( dot_step, -1.0 );
    sol->PlusAX( dot_step, (-1.0) * gamma * dt );

    dot_sol_alpha.PlusAX( dot_step, (-1.0) * alpha_m );
    sol_alpha.PlusAX( dot_step, (-1.0) * alpha_f * gamma * dt );


    // Assembly residual (& tangent if condition satisfied)
    if( nl_counter % nrenew_freq == 0 || nl_counter >= nrenew_threshold )
    {

      if(gbc-> UserLPM_Dirichlet_flag){
      // call GenBC to update Dirichlet and Neumann boundary conditions
        get_gbc_pressure_flow(curr_time,dt,alpha_f,sol_base,sol,dot_sol,&sol_alpha,infnbc_part,ebc_part,gbc,lassem_ptr, gassem_ptr,elements,quad_s,true);

      }

      gassem_ptr->Clear_KG();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(mat_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_tangent_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part, gbc );

#ifdef PETSC_USE_LOG
      PetscLogEventEnd(mat_assem_1_event,0,0,0,0);
#endif

      SYS_T::commPrint("  --- M updated");
      lsolver_ptr->SetOperator(gassem_ptr->K);
    }
    else
    {

      if(gbc-> UserLPM_Dirichlet_flag){
      // call GenBC to update Dirichlet and Neumann boundary conditions
        get_gbc_pressure_flow(curr_time,dt,alpha_f,sol_base,sol,dot_sol,&sol_alpha,infnbc_part,ebc_part,gbc,lassem_ptr, gassem_ptr,elements,quad_s,false);

      }

      gassem_ptr->Clear_G();

#ifdef PETSC_USE_LOG
      PetscLogEventBegin(vec_assem_1_event, 0,0,0,0);
#endif

      gassem_ptr->Assem_residual( &dot_sol_alpha, &sol_alpha, dot_sol, sol,
          curr_time, dt, alelem_ptr, lassem_ptr, elementv, elements,
          quad_v, quad_s, lien_ptr, anode_ptr,
          feanode_ptr, nbc_part, ebc_part, gbc );

#ifdef PETSC_USE_LOG
      PetscLogEventEnd(vec_assem_1_event,0,0,0,0);
#endif
    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);

    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);

    relative_error = residual_norm / initial_norm;

    if( relative_error >= nd_tol )
      SYS_T::print_fatal( "Error: nonlinear solver is diverging with error %e. Job killed.\n", relative_error);

  }while(nl_counter<nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol) conv_flag = true;
  else conv_flag = false;
}


void PNonlinear_NS_Solver::rescale_inflow_value( const double &stime,
    const ALocal_Inflow_NodalBC * const &infbc,
    const ICVFlowRate * const &flrate,
    const PDNSolution * const &sol_base,
    PDNSolution * const &sol ) const
{
  const int numnode = infbc -> get_Num_LD();

  const double val = flrate -> get_flow_rate( stime );
  double base_vals[3];
  int base_idx[3];

  for(int ii=0; ii<numnode; ++ii)
  {
    const int node_index = infbc -> get_LDN( ii );

    base_idx[0] = node_index * 4 + 1;
    base_idx[1] = node_index * 4 + 2;
    base_idx[2] = node_index * 4 + 3;

    VecGetValues(sol_base->solution, 3, base_idx, base_vals);

    VecSetValue(sol->solution, node_index*4+1, base_vals[0] * val, INSERT_VALUES);
    VecSetValue(sol->solution, node_index*4+2, base_vals[1] * val, INSERT_VALUES);
    VecSetValue(sol->solution, node_index*4+3, base_vals[2] * val, INSERT_VALUES);
  }

  VecAssemblyBegin(sol->solution); VecAssemblyEnd(sol->solution);
  sol->GhostUpdate();
}


void PNonlinear_NS_Solver::impose_inflow_value( const double &flow_rate,
    const ALocal_Inflow_NodalBC * const &infbc,
    const PDNSolution * const &sol_base,
    PDNSolution * const &sol ) const
{
  const int numnode = infbc -> get_Num_LD();

  const double val = flow_rate;

  double base_vals[3];
  int base_idx[3];

  for(int ii=0; ii<numnode; ++ii)
  {
    const int node_index = infbc -> get_LDN( ii );

    base_idx[0] = node_index * 4 + 1;
    base_idx[1] = node_index * 4 + 2;
    base_idx[2] = node_index * 4 + 3;

    VecGetValues(sol_base->solution, 3, base_idx, base_vals);

    VecSetValue(sol->solution, node_index*4+1, base_vals[0] * val, INSERT_VALUES);
    VecSetValue(sol->solution, node_index*4+2, base_vals[1] * val, INSERT_VALUES);
    VecSetValue(sol->solution, node_index*4+3, base_vals[2] * val, INSERT_VALUES);
  }

  VecAssemblyBegin(sol->solution); VecAssemblyEnd(sol->solution);
  sol->GhostUpdate();
}




void PNonlinear_NS_Solver::get_gbc_pressure_flow( const double &stime,
    const double &dt,
    const double &alpha_f,
    const PDNSolution * const &sol_base,
    PDNSolution * const &sol,
    PDNSolution * const &dot_sol,
    PDNSolution * const &sol_alpha,
    const ALocal_Inflow_NodalBC * const &infnbc,
    const ALocal_EBC * const &ebc_part,
    IGenBC * const &gbc,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s,
    const bool &compute_m_n_flag ) const
{

  int num_Neumann_faces=ebc_part->get_num_ebc();
  double * dot_flrate=new double[num_Neumann_faces];
  double * flrate=new double[num_Neumann_faces];
  double * P_n=new double[num_Neumann_faces];
  double * P_np1=new double[num_Neumann_faces];


  int num_Dirichlet_faces=gbc->get_num_Dirichlet_faces();


  double * inlet_avepre=new double[num_Dirichlet_faces];

  double * inlet_flrate_curr=new double[num_Dirichlet_faces];
  double * inlet_flrate_new=new double[num_Dirichlet_faces];

  double * inlet_flrate_alpha=new double[num_Dirichlet_faces];

  for(int ebc_id = 0; ebc_id < num_Neumann_faces; ++ebc_id)
  {
    // Calculate dot flow rate for face with ebc_id and MPI_Allreduce them
    // Here, dot_sol is the solution at time step n+1 (not n+alpha_f!)
    dot_flrate [ebc_id] = gassem_ptr->Assem_surface_flowrate( dot_sol, lassem_ptr,
        element_s, quad_s, ebc_part, ebc_id );

    // Calculate flow rate for face with ebc_id and MPI_Allreduce them
    // Here, sol is the solution at time step n+1 (not n+alpha_f!)
    flrate [ebc_id] = gassem_ptr->Assem_surface_flowrate( sol, lassem_ptr,
        element_s, quad_s, ebc_part, ebc_id );
     //   printf("outlet=%d flrate=%lf dot_flrate=%lf \n",ebc_id,flrate[ebc_id],dot_flrate[ebc_id]);

  }

  //currentlt 1 Dirichlet inlet only
  for(int ii=0;ii<num_Dirichlet_faces;++ii){
    inlet_avepre[ii]=gassem_ptr -> Assem_surface_ave_pressure(
      sol, lassem_ptr, element_s, quad_s, infnbc );
  }

  gbc->get_Q0(inlet_flrate_curr);
  const bool output_alldata_flag=false;
  gbc->get_P_Q(dot_flrate, flrate,inlet_avepre,P_np1,inlet_flrate_new,output_alldata_flag);


  for(int ii=0;ii<num_Neumann_faces;++ii){
     gbc->set_curr_P(ii,P_np1[ii]);
  }

  for(int ii=0;ii<num_Dirichlet_faces;++ii){
     gbc->set_curr_Q(ii,inlet_flrate_new[ii]);
  }



  for(int ii=0;ii<num_Dirichlet_faces;++ii){
   inlet_flrate_alpha[ii]=inlet_flrate_curr[ii]+(inlet_flrate_new[ii]-inlet_flrate_curr[ii])/dt*alpha_f;

   impose_inflow_value(inlet_flrate_new[ii],infnbc,sol_base,sol);
   impose_inflow_value(inlet_flrate_alpha[ii],infnbc,sol_base,sol_alpha);

  }


  if(compute_m_n_flag){

   double * m=new double[num_Neumann_faces];

   gbc->get_m(dot_flrate,flrate,inlet_avepre,m);

   for(int ii =0; ii<num_Neumann_faces;++ii){
     gbc->set_curr_m(ii,m[ii]);
     gbc->set_curr_n(ii,0.0);
   }

   delete [] m;m=nullptr;



  }



  /*
  //testing
    //  printf("time=%lf pressure=%lf \n",stime,inlet_face_avepre);


    // Get the (pressure) value on the outlet surface for traction evaluation
     gbc -> get_P0(P_n);
     gbc -> get_P(dot_flrate, flrate,P_np1 );

     printf("get_gbc_pressure_flow, p_np1=%lf %lf \n",P_np1[0],P_np1[1]);
    // Get the (potentially approximated) m := dP/dQ
     gbc -> get_m(dot_flrate, flrate,m_val );
     printf("get_gbc_pressure_flow, m_val=%lf %lf \n",m_val[0],m_val[1]);
    // Get the (potentially approximated) n := dP/d(dot_Q)
     gbc -> get_n( dot_flrate, flrate,n_val );
     printf("get_gbc_pressure_flow, n_val=%lf %lf \n",n_val[0],n_val[1]);
   */



  delete [] dot_flrate; dot_flrate=nullptr;
  delete [] flrate; flrate=nullptr;
  delete [] P_n; P_n=nullptr;
  delete [] P_np1; P_np1=nullptr;

  delete [] inlet_flrate_curr; inlet_flrate_curr=nullptr;
  delete [] inlet_flrate_new; inlet_flrate_new=nullptr;
  delete [] inlet_flrate_alpha; inlet_flrate_alpha=nullptr;
  delete [] inlet_avepre; inlet_avepre=nullptr;

}
// EOF
