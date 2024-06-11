#ifndef INTERFACE_PAIR_HPP
#define INTERFACE_PAIR_HPP

#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "VTK_Tools.hpp"
#include "HDF5_Tools.hpp"
#include "IIEN.hpp"
#include "Tet_Tools.hpp"
#include "Hex_Tools.hpp"

class Interface_pair
{
  public:
    // Constructor for type 0 interface-pair
    Interface_pair( const std::string &fixed_vtkfile,
                    const std::string &rotated_vtkfile,
                    const std::string &fixed_h5file,
                    const int &total_num_fixed_elem,
                    const int &total_num_fixed_pt,
                    const std::vector<double> &all_vol_ctrlPts,
                    const IIEN * const &VIEN,
                    const int &elemtype_in,
                    const std::vector<double> &intervals_in,
                    const int &direction_in );

    // Constructor for type 1 interface-pair
    Interface_pair( const std::string &fixed_vtkfile,
                    const std::string &rotated_vtkfile,
                    const std::string &fixed_h5file,
                    const int &total_num_fixed_elem,
                    const int &total_num_fixed_pt,
                    const std::vector<double> &all_vol_ctrlPts,
                    const IIEN * const &VIEN,
                    const int &elemtype_in,
                    const std::vector<double> &intervals_in,
                    const Vector_3 &centroid_in );

    virtual int get_num_fixed_ele() const
    {return  num_fixed_ele;}

    virtual int get_fixed_part_tag(const int &cell_index) const
    {return fixed_part_tag[cell_index];}

    virtual int get_fixed_faceID(const int &cell_index) const
    {return fixed_face_id[cell_index];}

    virtual std::vector<int> get_FL_vien() const
    {return fixed_layer_vien;}

    virtual std::vector<int> get_FLN_GID() const
    {return fixed_layer_global_node;}

    virtual std::vector<double> get_FLN_xyz() const
    {return fixed_layer_pt_xyz;}

    virtual std::vector<int> get_FIT() const
    {return fixed_interval_tag;}

    virtual std::vector<int> get_rotated_faceID() const
    {return rotated_face_id;}

    virtual std::vector<int> get_RL_vien() const
    {return rotated_layer_vien;}

    virtual std::vector<int> get_RLN_GID() const
    {return rotated_layer_global_node;}

    virtual std::vector<double> get_RLN_xyz() const
    {return rotated_layer_pt_xyz;}

    virtual std::vector<int> get_RIT() const
    {return rotated_interval_tag;}

    virtual ~Interface_pair(){};

  private:
    const int interface_type;
    // 0: Lofted along an axis
    // 1: Top or bottom

    int s_nLocBas;
    int v_nLocBas;

    const int T0_axial_direction;
    const Vector_3 T1_surface_centroid;

    int num_fixed_ele;
    std::vector<int> fixed_part_tag;
    std::vector<int> fixed_face_id;
    std::vector<int> fixed_layer_vien;
    std::vector<int> fixed_layer_global_node;
    std::vector<double> fixed_layer_pt_xyz;
    std::vector<int> fixed_interval_tag;

    int num_rotated_ele;
    std::vector<int> rotated_face_id;
    std::vector<int> rotated_layer_vien;
    std::vector<int> rotated_layer_global_node;
    std::vector<double> rotated_layer_pt_xyz;
    std::vector<int> rotated_interval_tag;

    virtual void Initialize(const std::string &fixed_vtkfile,
      const std::string &rotated_vtkfile,
      const std::string &fixed_h5file,
      const int &total_num_fixed_elem,
      const int &total_num_fixed_pt,
      const std::vector<double> &all_vol_ctrlPts,
      const IIEN * const &VIEN,
      const int &elemtype_in,
      const std::vector<double> &intervals_in);

    virtual void Check_interval(const std::vector<double> &intervals);

    virtual void Tag(const std::vector<double> &intervals,
      const int &num_fixed_ele,
      const std::vector<double> &fixed_sur_pt_xyz,
      const std::vector<int> &fixed_sur_ien,
      const int &num_rotated_ele,
      const std::vector<double> &rotated_sur_pt_xyz,
      const std::vector<int> &rotated_sur_ien);

    virtual int Group(const Vector_3 &ele_centroid, const std::vector<double> &intervals);

    Interface_pair() = delete;

};

#endif