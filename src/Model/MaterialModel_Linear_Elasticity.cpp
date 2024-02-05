#include "MaterialModel_Linear_Elasticity.hpp"

MaterialModel_Linear_Elasticity::MaterialModel_Linear_Elasticity(
    const double &in_modulus_E, const double &in_nu )
: modulus_E( in_modulus_E ), nu( in_nu ),
  lambda( in_nu * in_modulus_E / ((1.0 + in_nu) * (1.0 - 2.0 * in_nu)) ),
  mu( 0.5 * in_modulus_E / (1.0 + in_nu) )
{
}

MaterialModel_Linear_Elasticity::MaterialModel_Linear_Elasticity(
    const char * const &fname)
{
  hid_t h5file = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * h5r = new HDF5_Reader( h5file );

  SYS_T::print_fatal_if( h5r->read_string("/", "model_name") != get_model_name(),
     "Error: MaterialModel_Linear_Elasticity constructor does not match h5 file.\n" );

  modulus_E = h5r -> read_doubleScalar("/", "modulus_E");
  nu        = h5r -> read_doubleScalar("/", "nu");
  lambda    = h5r -> read_doubleScalar("/", "lambda");
  mu        = h5r -> read_doubleScalar("/", "mu");
  
  delete h5r; H5Fclose(h5file);
}

MaterialModel_Linear_Elasticity::~MaterialModel_Linear_Elasticity()
{}

void MaterialModel_Linear_Elasticity::print_info() const
{
  SYS_T::commPrint("\t  MaterialModel_Linear_Elasticity: \n");
  SYS_T::commPrint("\t  Young's Modulus E     = %e \n", modulus_E);
  SYS_T::commPrint("\t  Possion's ratio nu    = %e \n", nu);
  SYS_T::commPrint("\t  Lame parameter lambda = %e \n", lambda);
  SYS_T::commPrint("\t  Lame parameter mu     = %e \n", mu);
}

void MaterialModel_Linear_Elasticity::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * h5w = new HDF5_Writer(file_id);

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar("modulus_E", modulus_E);
    h5w -> write_doubleScalar("nu", nu);
    h5w -> write_doubleScalar("lambda", lambda);
    h5w -> write_doubleScalar("mu", mu);

    delete h5w; H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

Tensor2_3D MaterialModel_Linear_Elasticity::get_Cauchy_stress( const Tensor2_3D &F ) const
{
  const double l2mu = lambda + 2.0 * mu;

  Tensor2_3D stress;
  stress(0) = l2mu * F(0) + lambda * (F(4) + F(8));
  stress(1) = mu * ( F(1) + F(3) );
  stress(2) = mu * ( F(2) + F(6) );
  stress(3) = mu * ( F(1) + F(3) );
  stress(4) = l2mu * F(4) + lambda * (F(0) + F(8));
  stress(5) = mu * ( F(5) + F(7) );
  stress(6) = mu * ( F(2) + F(6) );
  stress(7) = mu * ( F(5) + F(7) );
  stress(8) = l2mu * F(8) + lambda * (F(0) + F(4));

  return stress;
}

// EOF