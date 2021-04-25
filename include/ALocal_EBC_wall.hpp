#ifndef ALOCAL_EBC_WALL_HPP
#define ALOCAL_EBC_WALL_HPP
// ==================================================================
// ALocal_EBC_wall.hpp
//
// Analysis use: local subdomain's elemental boundary condition. This
// is a derived class from ALocal_EBC to add the wall thickness and
// young's modulus, as well as the fluid density used to compute the
// young's modulus.
//
// Author: Ju Liu
// Date: Aug. 10 2020
// ==================================================================
#include "ALocal_EBC.hpp"
#include "HDF5_Writer.hpp"

class ALocal_EBC_wall : public ALocal_EBC
{
  public:
    ALocal_EBC_wall( const std::string &fileBaseName,
        const int &in_cpu_rank, const IQuadPts * const &quad,
        const std::string &gname = "/ebc",
        const bool &prestress_flag = false );

    virtual ~ALocal_EBC_wall();

    virtual void print_info() const;

    // ee is the local element index, ranges in [ 0, num_local_cell[0] )
    // elem_thickness is the output with length cell_nLocBas[0], which
    // represents the thickness for each local node
    virtual void get_thickness(const int &ee, double * const &elem_thickness) const;

    // elem_youngsmod is the output with length cell_nLocBas[0], which
    // represents the Young's modulus at each local node
    virtual void get_youngsmod(const int &ee, double * const &elem_youngsmod) const;

    // elem_quaprestress is the output with length 6 x num_local_cell[0] x
    // num_face_nqp, which represents the prestress values at each quadrature
    // points. The data is arranged the following way,
    // elem_quaprestress[6*ii+jj] is the stress tensor jj-th component in Voigt
    // notation at the ii-th quadrature point.
    virtual void get_prestress(const int &ee, double * const &elem_quaprestress) const;

    // elem_quaprestress is the input with length 6 x num_face_nqp.
    virtual void set_prestress(const int &ee, const double * const &elem_quaprestress);

    virtual double get_fluid_density() const {return fluid_density;}

    virtual void write_prestress_hdf5( const char * FileName ) const;

  protected:
    // The rank or id of the subdomain
    const int cpu_rank;

    // num quadrature points for surface element
    const int face_nqp;

    // flag for whether the wall prestress is being solved for and updated 
    const bool solve_prestress;

    // Fluid density used to generate the youngs modulus for arteries at
    // the preprocessing stage.
    double fluid_density;

    // If this partition owns any part of the wall, the thickness and
    // youngsmod vectors are each of length num_local_node[0].
    // Otherwise, these vectors are of length 0.
    std::vector<double> thickness, youngsmod;

    // Prestress components (6 in Voigt notation) per quadrature point
    // Length 6 x face_nqp x num_local_cell[0]
    std::vector<double> qua_prestress;
};

#endif
