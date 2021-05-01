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
  
  Qi0 = 0.0;
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
    const IGenBC * const &gbc,
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

  // ------------------------------------------------- 
  // Update the inflow boundary values
  //rescale_inflow_value(curr_time+dt, infnbc_part, flr_ptr, sol_base, sol);
  //rescale_inflow_value(curr_time+alpha_f*dt, infnbc_part, flr_ptr, sol_base, &sol_alpha);
  // ------------------------------------------------- 
  
  
  set_Dirichlet_inflow( curr_time,dt,gamma,alpha_f,alpha_m,sol_base,pre_sol,pre_dot_sol,
      sol,dot_sol,&sol_alpha,&dot_sol_alpha,infnbc_part,lassem_ptr, gassem_ptr,elements,quad_s );
  
  

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
    
    set_Dirichlet_inflow( curr_time,dt,gamma,alpha_f,alpha_m,sol_base,pre_sol,pre_dot_sol,
        sol,dot_sol,&sol_alpha,&dot_sol_alpha,infnbc_part, lassem_ptr, gassem_ptr,elements,quad_s );    

    // Assembly residual (& tangent if condition satisfied) 
    if( nl_counter % nrenew_freq == 0 || nl_counter >= nrenew_threshold )
    {
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
    
    SYS_T::print_fatal_if( residual_norm != residual_norm, "Error: nonlinear solver residual norm is NaN. Job killed.\n" );
    
    SYS_T::commPrint("  --- nl_res: %e \n", residual_norm);
    
    if( initial_norm < 1e-7 ) initial_norm = 1e-7;
    relative_error = residual_norm / initial_norm;

    SYS_T::print_fatal_if( relative_error >= nd_tol, "Error: nonlinear solver is diverging with error %e. Job killed.\n", relative_error);

  }while(nl_counter<nmaxits && relative_error > nr_tol && residual_norm > na_tol);

  Print_convergence_info(nl_counter, relative_error, residual_norm);
  
  update_Qi0();
  
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


void PNonlinear_NS_Solver::set_Dirichlet_inflow( const double &stime,
    const double &dt,
    const double &gamma,
    const double &alpha_f,
    const double &alpha_m,
    const PDNSolution * const &sol_base,
    const PDNSolution * const &pre_sol,
    const PDNSolution * const &pre_dot_sol,
    PDNSolution * const &sol,
    PDNSolution * const &dot_sol,
    PDNSolution * const &sol_alpha,
    PDNSolution * const &dot_sol_alpha,
    const ALocal_Inflow_NodalBC * const &infnbc,
    IPLocAssem * const &lassem_ptr,
    IPGAssem * const &gassem_ptr,
    FEAElement * const &element_s,
    const IQuadPts * const &quad_s ) const
{
  //currentlt 1 Dirichlet inlet only
  const double inlet_avepre = gassem_ptr -> Assem_surface_ave_pressure(
      sol, lassem_ptr, element_s, quad_s, infnbc );
  
  const double inlet_avepre_prev = gassem_ptr -> Assem_surface_ave_pressure(
      pre_sol, lassem_ptr, element_s, quad_s, infnbc );
  // use last 0D solution as initial Q for the next 0D integration
  const double inlet_flrate_prev = Qi0;
  
  const double inlet_flrate = get_inductor_flow( inlet_avepre, inlet_avepre_prev, inlet_flrate_prev, stime, dt ); 

  // if(SYS_T::get_MPI_rank() == 0)
  //   std::cout << "t = " << stime << " P = " << inlet_avepre << " Q = " << inlet_flrate << "\n";

  const double inlet_flrate_alpha = inlet_flrate_prev + ( inlet_flrate - inlet_flrate_prev ) * alpha_f;

  impose_inflow_value_correct_ac( inlet_flrate,inlet_flrate_alpha,dt,gamma,alpha_f,alpha_m,infnbc,
      sol_base, pre_sol,pre_dot_sol,sol,dot_sol,sol_alpha,dot_sol_alpha );
}

double PNonlinear_NS_Solver::get_inductor_flow( const double &inlet_press, 
    const double &inlet_press_prev,
    const double &inlet_flow_prev,
    const double &curr,
    const double &dt ) const
{
  const double fac13 = 1.0 / 3.0;
  const double fac23 = 2.0 / 3.0;
  const double fac18 = 1.0 / 8.0;
  const double fac38 = 3.0 / 8.0;
  const int N = 1000;
  const double h = dt / static_cast<double>(N);

  double Qi_m = inlet_flow_prev; 
  
  for(int mm=0; mm<N; ++mm)
  {
    const double P_m = inlet_press_prev + static_cast<double>(mm) * ( inlet_press - inlet_press_prev ) / static_cast<double>(N);
    
    const double P_mp1 = inlet_press_prev + static_cast<double>(mm+1) * ( inlet_press - inlet_press_prev ) / static_cast<double>(N);

    const double K1 = F( P_m, curr );

    const double K2 = F( fac23 * P_m + fac13 * P_mp1, curr );

    const double K3 = F( fac13 * P_m + fac23 * P_mp1, curr );
  
    const double K4 = F( P_mp1, curr );

    Qi_m = Qi_m + fac18 * K1 * h + fac38 * K2 * h + fac38 * K3 * h + fac18 * K4 * h; 
  }

  // if(SYS_T::get_MPI_rank() == 0)
  //   std::cout << "\ninlet_flow_prev =" << inlet_flow_prev << " Qi_m = " << Qi_m << " F =" << F(inlet_press, curr) << "\n";
 
  record_0D_Q(Qi_m);
  
  return Qi_m;
  //return -Qi_m + F( inlet_press , curr) * dt;
}    

double PNonlinear_NS_Solver::F( const double &Pi, const double &curr ) const
{
  const double P_head = 5000.0;
  const double LL = 2.0;
  return ( P_head - Pi ) / LL;
}

void PNonlinear_NS_Solver::impose_inflow_value_correct_ac( const double &flow_rate,
    const double &flow_rate_alpha,
    const double &dt,
    const double &gamma,
    const double &alpha_f,
    const double &alpha_m,
    const ALocal_Inflow_NodalBC * const &infbc,
    const PDNSolution * const &sol_base,
    const PDNSolution * const &pre_sol,
    const PDNSolution * const &pre_dot_sol,
    PDNSolution * const &sol,
    PDNSolution * const &dot_sol,
    PDNSolution * const &sol_alpha,
    PDNSolution * const &dot_sol_alpha ) const
{
  const int numnode = infbc -> get_Num_LD();

  double base_vals[3];
  int base_idx[3];

  double pre_vals[3];
  double curr_vals[3];
  double pre_dot_vals[3];
  double dot_vals[3];

  for(int ii=0; ii<numnode; ++ii)
  {
    const int node_index = infbc -> get_LDN( ii );

    base_idx[0] = node_index * 4 + 1;
    base_idx[1] = node_index * 4 + 2;
    base_idx[2] = node_index * 4 + 3;

    VecGetValues(sol->solution, 3, base_idx, curr_vals);
    VecGetValues(sol_base->solution, 3, base_idx, base_vals);

    curr_vals[0]=base_vals[0]*flow_rate;
    curr_vals[1]=base_vals[1]*flow_rate;
    curr_vals[2]=base_vals[2]*flow_rate;

    VecSetValue(sol->solution, node_index*4+1, curr_vals[0], INSERT_VALUES);
    VecSetValue(sol->solution, node_index*4+2, curr_vals[1], INSERT_VALUES);
    VecSetValue(sol->solution, node_index*4+3, curr_vals[2], INSERT_VALUES);


    //make correction for acceleration values on Dirichal nodes
    // Ydot_n+1=(gamma-1)/gamma*Ydot_n+(Y_n+1 - Y_n)/(gamma*dt)

    VecGetValues(pre_sol->solution,3,base_idx,pre_vals);
    VecGetValues(pre_dot_sol->solution,3,base_idx,pre_dot_vals);
    VecGetValues(dot_sol->solution,3,base_idx,dot_vals);

    dot_vals[0]=(gamma-1.0)/gamma*pre_dot_vals[0]+(curr_vals[0]-pre_vals[0])/(gamma*dt);
    dot_vals[1]=(gamma-1.0)/gamma*pre_dot_vals[1]+(curr_vals[1]-pre_vals[1])/(gamma*dt);
    dot_vals[2]=(gamma-1.0)/gamma*pre_dot_vals[2]+(curr_vals[2]-pre_vals[2])/(gamma*dt);

    VecSetValue(dot_sol->solution, node_index*4+1, dot_vals[0], INSERT_VALUES);
    VecSetValue(dot_sol->solution, node_index*4+2, dot_vals[1], INSERT_VALUES);
    VecSetValue(dot_sol->solution, node_index*4+3, dot_vals[2], INSERT_VALUES);

    curr_vals[0]=(1.0-alpha_f)*pre_vals[0]+alpha_f*curr_vals[0];
    curr_vals[1]=(1.0-alpha_f)*pre_vals[1]+alpha_f*curr_vals[1];
    curr_vals[2]=(1.0-alpha_f)*pre_vals[2]+alpha_f*curr_vals[2];
    
    VecSetValue(sol_alpha->solution, node_index*4+1, curr_vals[0], INSERT_VALUES);
    VecSetValue(sol_alpha->solution, node_index*4+2, curr_vals[1], INSERT_VALUES);
    VecSetValue(sol_alpha->solution, node_index*4+3, curr_vals[2], INSERT_VALUES);

    dot_vals[0]=(1.0-alpha_m)*pre_dot_vals[0]+alpha_m*dot_vals[0];
    dot_vals[1]=(1.0-alpha_m)*pre_dot_vals[1]+alpha_m*dot_vals[1];
    dot_vals[2]=(1.0-alpha_m)*pre_dot_vals[2]+alpha_m*dot_vals[2];
    
    VecSetValue(dot_sol_alpha->solution, node_index*4+1, dot_vals[0], INSERT_VALUES);
    VecSetValue(dot_sol_alpha->solution, node_index*4+2, dot_vals[1], INSERT_VALUES);
    VecSetValue(dot_sol_alpha->solution, node_index*4+3, dot_vals[2], INSERT_VALUES);
  }

  VecAssemblyBegin(sol->solution); VecAssemblyEnd(sol->solution);
  sol->GhostUpdate();

  VecAssemblyBegin(dot_sol->solution); VecAssemblyEnd(dot_sol->solution);
  dot_sol->GhostUpdate();
  
  VecAssemblyBegin(sol_alpha->solution); VecAssemblyEnd(sol_alpha->solution);
  sol_alpha->GhostUpdate();
  
  VecAssemblyBegin(dot_sol_alpha->solution); VecAssemblyEnd(dot_sol_alpha->solution);
  dot_sol_alpha->GhostUpdate();
}

void PNonlinear_NS_Solver::record_0D_Q( const double & Qim) const
{
  Qim_final = Qim;
}

void PNonlinear_NS_Solver::update_Qi0() const
{
  Qi0 = Qim_final;
}

// EOF
