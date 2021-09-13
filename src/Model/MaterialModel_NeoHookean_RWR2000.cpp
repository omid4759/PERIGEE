#include "MaterialModel_NeoHookean_RWR2000.hpp"

MaterialModel_NeoHookean_RWR2000::MaterialModel_NeoHookean_RWR2000( 
    const double &in_E, const double &in_nu )
: E( in_E ), nu( in_nu ), lambda( nu * E / ((1+nu) * (1-2.0*nu)) ),
  mu( E/(2.0+2.0*nu) ), kappa( lambda + 2.0 * mu / 3.0 )
{}


MaterialModel_NeoHookean_RWR2000::~MaterialModel_NeoHookean_RWR2000()
{}


void MaterialModel_NeoHookean_RWR2000::print_info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  MaterialModel_NeoHookean_RWR2000: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  Young's Modulus E  = %e \n", E);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Possion's ratio nu = %e \n", nu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Lame coeff lambda  = %e \n", lambda);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Shear modulus mu   = %e \n", mu);
  PetscPrintf(PETSC_COMM_WORLD, "\t  Bulk modulus kappa = %e \n", kappa);
}

void MaterialModel_NeoHookean_RWR2000::write_hdf5( const char * const &fname ) const
{
  if( SYS_T::get_MPI_rank() == 0 )
  {
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_Writer * h5w = new HDF5_Writer(file_id);

    h5w -> write_string("model_name", get_model_name());
    h5w -> write_doubleScalar("E", E);
    h5w -> write_doubleScalar("nu", nu);
    h5w -> write_doubleScalar("lambda", lambda);
    h5w -> write_doubleScalar("mu", mu);
    h5w -> write_doubleScalar("kappa", kappa);

    delete h5w; H5Fclose(file_id);
  }

  MPI_Barrier(PETSC_COMM_WORLD);
}

void MaterialModel_NeoHookean_RWR2000::get_PK( 
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S )
{
  Cinv.MatMultTransposeLeft(F);
  Cinv.inverse();
  const double ln_detF = std::log( F.det() );
  S.gen_id(); S.scale( mu );
  S.AXPY( kappa * ln_detF - mu , Cinv );
  P.MatMult( F, S );
}


void MaterialModel_NeoHookean_RWR2000::get_PK_Stiffness( 
    const Matrix_3x3 &F, Matrix_3x3 &P, Matrix_3x3 &S, Tensor4_3D &CC)
{
  Cinv.MatMultTransposeLeft(F);
  Cinv.inverse();
  const double ln_detF = std::log( F.det() );
  
  S.gen_id(); S.scale( mu );
  S.AXPY( kappa * ln_detF - mu , Cinv );
  P.MatMult( F, S );

  CC.gen_zero();
  const double val1 = 2.0 * ( mu - kappa * ln_detF );
  CC.add_OutProduct(kappa, Cinv, Cinv);
  CC.add_SymmProduct(val1, Cinv, Cinv);
}


double MaterialModel_NeoHookean_RWR2000::get_strain_energy( const Matrix_3x3 &F )
{
  Cinv.MatMultTransposeLeft(F);
  const double trC = Cinv.tr();
  const double lnJ = std::log( F.det() );
  return 0.5 * mu * (trC - 3.0) - mu * lnJ + 0.5 * kappa * lnJ * lnJ;
}

// EOF
