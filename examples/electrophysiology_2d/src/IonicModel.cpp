#include "IonicModel.hpp"

IonicModel::IonicModel()
  :ap_1{100}, ap_2{80}, ap_3{12.9}, m1{0.2}, m2{0.3},alpha{0.01},
   gamma{0.002}, b{0.15}, c{8}, d_iso{0.1}, d_ani{0.0}, tol{1e-8},
   chi{140}, C_m{1}
{
  SYS_T::commPrint("IonicModel constructor. \n");
};

IonicModel::~IonicModel()
{
  SYS_T::commPrint("IonicModel destructor. \n");
};

void IonicModel::print_info () const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  Aliev Panfilov EP: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_1 = %e \n", ap_1);
  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_2 = %e \n", ap_2);
  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_3 = %e \n", ap_3);
};

double IonicModel::get_diso() const
{
  return d_iso;
};

double IonicModel::get_dani() const
{
  return d_ani;
};

double IonicModel::get_chi() const
{
  return chi;
};

double IonicModel::get_C_m() const
{
  return C_m;
};

void IonicModel::material_routine(const double &r_old_in,
				  const double &dt_in,
				  const double &Phi_in,
				  double &f_Phi,
				  double &dP_fP,
				  double &r_new) const
{
//  double dP_dr=0;
//  //non dimensionalize 
//  const double Phi_nd { (Phi_in+ap_2)/ap_1};
//  const double dt_nd { dt_in/ap_3};
//
//  //find r_new and dP_dr from Newton Raphson iterations
//  double Rr, dr_Rr, dp_Rr;
//  unsigned int i_counter=0;
//  r_new=r_old_in; //initialize
//  Rr= r_new-r_old_in-((gamma+(m1*r_new)/(m2+Phi_nd))
//		   *(-r_new-c*Phi_nd*(Phi_nd-b-1)))*dt_nd;
//
//  while (std::abs(Rr)>tol){
//    ++i_counter;
//    if (i_counter>=10) {
//      std::cout << "no convergence at local newton iteration: "
//		<< i_counter << std::endl;
//      throw std::runtime_error("No convergence.");
//    }
//    dr_Rr=1+(gamma+m1/(m2+Phi_nd)*(2*r_new+c*Phi_nd*(Phi_nd-b-1)))*dt_nd;
//    r_new=r_new-Rr/dr_Rr;//update r_new
//    Rr= r_new-r_old_in-((gamma+(m1*r_new)/(m2+Phi_nd))
//		     *(-r_new-c*Phi_nd*(Phi_nd-b-1)))*dt_nd;
//  }
//  
//  dr_Rr=1+(gamma+m1/(m2+Phi_nd)*(2*r_new+c*Phi_nd*(Phi_nd-b-1)))*dt_nd;
//  dp_Rr=((gamma+m1*r_new/(m2+Phi_nd))*c*(2*Phi_nd-b-1)-m1*r_new/pow((m2+Phi_nd),2)
//	 *(r_new+c*Phi_nd*(Phi_nd-b-1)))*dt_nd;
//  dP_dr=-dp_Rr/dr_Rr;
//
//  //r_new and dP_dr found now. calculate rest.
//  f_Phi = c*Phi_nd*(Phi_nd- alpha)*(1-Phi_nd)-r_new*Phi_nd;
//  dP_fP = c*(-3*pow(Phi_nd,2)
//	     +2*(1+alpha)*Phi_nd-alpha)-r_new-Phi_nd*dP_dr;

  //redimensionalize
  // f_Phi=ap_1*f_Phi/ap_3;
  //dP_fP=dP_fP/ap_3;
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //for manufactured solution:
  r_new=dt_in;
  f_Phi=10;
  dP_fP=0;
  
}


// EOF
