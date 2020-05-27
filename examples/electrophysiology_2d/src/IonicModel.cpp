#include "IonicModel.hpp"

IonicModel::IonicModel()
{
  SYS_T::commPrint("IonicModel constructor. \n");
};

IonicModel::~IonicModel()
{
  SYS_T::commPrint("IonicModel destructor. \n");
};

void IonicModel::print_info () const
{
  SYS_T::commPrint("IonicModel Info: NA. \n");
};

void IonicModel::material_routine(const double &r_old,
				  const double  &dt,
				  double &r_new_out) const
{
  r_new_out= r_old+dt;  // fill-in with the actual ioni model
}


// EOF
