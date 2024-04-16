#include "FEAElement_Triangle3_3D_der0.hpp"

FEAElement_Triangle3_3D_der0::FEAElement_Triangle3_3D_der0( 
    const int &in_nqua ) : numQuapts( in_nqua )
{
  R = new double [ 3 * numQuapts ];
}

FEAElement_Triangle3_3D_der0::~FEAElement_Triangle3_3D_der0()
{
  delete [] R; R = nullptr;
}

void FEAElement_Triangle3_3D_der0::print_info() const
{
  SYS_T::commPrint("Triangle3_3D_der0: ");
  SYS_T::commPrint("3-node triangle element with no derivative evaluated. \n ");
  PetscPrintf(PETSC_COMM_WORLD, "elemType: %d. \n", get_Type());
  SYS_T::commPrint("Note: This element is designed for natural BC integrals. \n ");
}

double FEAElement_Triangle3_3D_der0::get_memory_usage() const
{
  double double_size = 3 * numQuapts + 10.0;
  double int_size = 2;
  return double_size * 8.0 + int_size * 4.0;
}

void FEAElement_Triangle3_3D_der0::buildBasis( const IQuadPts * const &quad,
    const double * const &ctrl_x,
    const double * const &ctrl_y,
    const double * const &ctrl_z )
{
  for( int qua = 0; qua < numQuapts; ++qua )
  {
    const double qua_r = quad -> get_qp(qua, 0);
    const double qua_s = quad -> get_qp(qua, 1);
    R[qua*3 + 0] = 1.0 - qua_r - qua_s;
    R[qua*3 + 1] = qua_r;
    R[qua*3 + 2] = qua_s;
  }

  dx_dr = Vector_3( - ctrl_x[0] + ctrl_x[1],
                    - ctrl_y[0] + ctrl_y[1],
                    - ctrl_z[0] + ctrl_z[1]);

  dx_ds = Vector_3( - ctrl_x[0] + ctrl_x[2],
                    - ctrl_y[0] + ctrl_y[2],
                    - ctrl_z[0] + ctrl_z[2]);

  // vec(un) = vec(dx_dr) x vec(dx_ds)
  un = Vec3::cross_product( dx_dr, dx_ds );

  // area = || vec(un) ||
  detJac = un.normalize();
}

void FEAElement_Triangle3_3D_der0::get_R( const int &quaindex, 
    double * const &basis ) const
{
  ASSERT(quaindex>=0 && quaindex < numQuapts, "FEAElement_Triangle3_3D_der0::get_R function error.\n" );
  const int offset = quaindex * 3;
  basis[0] = R[offset];
  basis[1] = R[offset+1];
  basis[2] = R[offset+2];
}

std::vector<double> FEAElement_Triangle3_3D_der0::get_R( const int &quaindex ) const
{
  ASSERT(quaindex>=0 && quaindex < numQuapts, "FEAElement_Triangle3_3D_der0::get_R function error.\n" );
  const int offset = quaindex * 3;
  return { R[offset], R[offset+1], R[offset+2] };
}

Vector_3 FEAElement_Triangle3_3D_der0::get_2d_normal_out( const int &quaindex,
    double &area ) const
{
  ASSERT(quaindex>=0 && quaindex < numQuapts, "FEAElement_Triangle3_3D_der0::get_2d_normal_out function error.\n" );
  area = detJac;
  return un;
}

Vector_3 FEAElement_Triangle3_3D_der0::get_normal_out( const int &quaindex,
    const Vector_3 &sur_pt, const Vector_3 &int_pt, double &area ) const
{
  // Construct a vector from the interior point to the triangle first node
  const Vector_3 mm = sur_pt - int_pt;

  // Dot product of the defined vector with the calculated normal vector
  const double mdotn = Vec3::dot_product( mm, un );

  SYS_T::print_fatal_if( std::abs(mdotn) < 1.0e-10, "Warning: FEAElement_Triangle3_3D_der0::get_normal_out, the element might be ill-shaped.\n");

  area = detJac;

  // If dot product is negative, adjust the normal vector
  if(mdotn < 0) return (-1.0)*un;
  else return un;
}

bool FEAElement_Triangle3_3D_der0::search_closed_point(
  IQuadPts * const &closest_point, const Vector_3 &target_xyz,
  const double * const &ctrl_x, const double * const &ctrl_y, const double * const &ctrl_z )
{
  // initial value
  buildBasis(closest_point, ctrl_x, ctrl_y, ctrl_z);

  std::vector<double> basis = get_R(0);

  Vector_3 point_xyz(0.0, 0.0, 0.0);
  for(int ii=0; ii<3; ++ii)
  {
    point_xyz(0) += basis[ii] * ctrl_x[ii];
    point_xyz(1) += basis[ii] * ctrl_y[ii];
    point_xyz(2) += basis[ii] * ctrl_z[ii];
  }

  // initial distance
  const double init_dist = (point_xyz - target_xyz).norm2();
  SYS_T::commPrint("init_dist: %e\n", init_dist);
  if (init_dist < 1e-9) return true;

  double curr_dist = init_dist;

  int iter_counter = 0;

  while(iter_counter < 100 && curr_dist > 1e-9)
  {
    const double Resr = 2 * (dx_dr(0) * (point_xyz(0) - target_xyz(0))
                         +   dx_dr(1) * (point_xyz(1) - target_xyz(1))
                         +   dx_dr(2) * (point_xyz(2) - target_xyz(2)));
  
    const double Ress = 2 * (dx_ds(0) * (point_xyz(0) - target_xyz(0))
                          +  dx_ds(1) * (point_xyz(1) - target_xyz(1))
                          +  dx_ds(2) * (point_xyz(2) - target_xyz(2)));

    const double dResr_dr = 2 * (dx_dr(0) * dx_dr(0) + dx_dr(1) * dx_dr(1) + dx_dr(2) * dx_dr(2));

    const double dRess_ds = 2 * (dx_ds(0) * dx_ds(0) + dx_ds(1) * dx_ds(1) + dx_ds(2) * dx_ds(2));

    const double dResr_ds = 2 * (dx_dr(0) * dx_ds(0) + dx_dr(1) * dx_ds(1) + dx_dr(2) * dx_ds(2));

    const double dRess_dr = dResr_ds;

    std::array<double, 4> d_mat = {dResr_dr, dResr_ds, dRess_dr, dRess_ds};
    MATH_T::Matrix_Dense<2> D_mat(d_mat);

    std::array<double, 2> Res_vec = {-Resr, -Ress};
    std::array<double, 2> dxi = D_mat.LU_solve(Res_vec);

    double xi_r = closest_point->get_qp(0, 0);
    double xi_s = closest_point->get_qp(0, 1);

    xi_r += dxi[0];
    xi_s += dxi[1];

    double xi_t = 1 - xi_r - xi_s;

    // Update the closest_point
    closest_point->set_qp(0, 0, xi_r);
    closest_point->set_qp(0, 1, xi_s);
    closest_point->set_qp(0, 2, xi_t);

    buildBasis(closest_point, ctrl_x, ctrl_y, ctrl_z);

    std::vector<double> basis = get_R(0);

    Vector_3 point_xyz(0.0, 0.0, 0.0);
    for(int ii=0; ii<3; ++ii)
    {
      point_xyz(0) += basis[ii] * ctrl_x[ii];
      point_xyz(1) += basis[ii] * ctrl_y[ii];
      point_xyz(2) += basis[ii] * ctrl_z[ii];
    }

    curr_dist = (point_xyz - target_xyz).norm2();
    SYS_T::commPrint("iter: %d, curr_dist: %e\n", iter_counter, curr_dist);

    iter_counter += 1;
  }

  if(curr_dist > 1e-9)
    return false;
   else
    return true;
}

// EOF
