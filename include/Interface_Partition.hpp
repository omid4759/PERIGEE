#ifndef INTERFACE_PARTITION_HPP
#define INTERFACE_PARTITION_HPP

#include "IPart.hpp"
#include "Map_Node_Index.hpp"
#include "HDF5_Writer.hpp"
#include "Interface_pair.hpp"

class Interface_Partition
{
  public:
    Interface_Partition( const IPart * const &part,
        const Map_Node_Index * const &mnindex,
        const std::vector<Interface_pair> &interfaces );

    virtual ~Interface_Partition(){};

    // write the data to hdf5 file in folder /interface
    virtual void write_hdf5( const std::string &FileName ) const;

  private:
    const int cpu_rank;

    const int num_pair;

    std::vector<int> num_tag;

    // the number of the local elements of the fixed interfaces
    std::vector<int> num_fixed_part_ele;

    // stores the face id of the volume element of the fixed interface in this part
    std::vector<std::vector<int>> fixed_ele_face_id;

    std::vector<std::vector<int>> fixed_layer_ien;
    std::vector<std::vector<int>> fixed_layer_global_node;
    std::vector<std::vector<double>> fixed_layer_pt_xyz;
    std::vector<std::vector<int>> fixed_interval_tag;

    // stores the face id of all the volume element of the rotated layer
    std::vector<std::vector<std::vector<int>>> rotated_ele_face_id;

    // stores ien of all the volume element of the rotated layer
    std::vector<std::vector<std::vector<int>>> rotated_layer_ien;

    // stores the rotated layer nodes indices, converted by get_old2new()
    std::vector<std::vector<int>> rotated_layer_global_node;

    // stores the rotated layer nodes' coordinates
    std::vector<std::vector<double>> rotated_layer_pt_xyz;

    std::vector<std::vector<int>> rotated_interval_tag;
};

#endif