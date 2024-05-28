#ifndef INTERFACE_ELEMENT_GROUP_HPP
#define INTERFACE_ELEMENT_GROUP_HPP

#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "VTK_Tools.hpp"

class Interface_Element_Group
{
  public:
    Interface_Element_Group(const std::string &vtkfilename, const int &nLocBas, const int &direction_in, const std::vector<double> &intervals_in);

    Interface_Element_Group(const std::string &vtkfilename, const int &nLocBas, const Vector_3 &centroid_in, const std::vector<double> &intervals_in);

    ~Interface_Element_Group() = default;

    std::vector<int> get_tag()
    {return interval_tag;}


  private:
    // Type:
    // 0 - intervals lay along an axis
    // 1 - intervals are concentric circles with increasing radii
    int groups_type;

    int axis_direction;

    Vector_3 centroid;

    std::vector<double> intervals;

    std::vector<int> interval_tag;

    void Tag(const std::string &vtkfilename, const int &nLocBas);

    int Group(const Vector_3 &ele_centroid);

    Interface_Element_Group() = delete;
};

#endif