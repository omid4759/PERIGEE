#include "Interface_Element_Group.hpp"

Interface_Element_Group::Interface_Element_Group(
  const std::string &vtkfilename, const int &nLocBas, const int &direction_in, const std::vector<double> &intervals_in) :
  groups_type {0}, axis_direction {direction_in}, centroid {Vector_3(0, 0, 0)},
  intervals {intervals_in}, interval_tag {std::vector<int> {}}
{
  // Check the intervals
  SYS_T::print_fatal_if(VEC_T::get_size(intervals) < 2,
    "Error, Interface_Element_Group: illegal 'interval' vector. (1)\n");

  const int n_interval = VEC_T::get_size(intervals) - 1;
  
  for(int ii{0}; ii<n_interval; ++ii)
  {
    SYS_T::print_fatal_if(intervals[ii] >= intervals[ii + 1],
      "Error, Interface_Element_Group: illegal 'interval' vector. (2)\n");
  }
  
  Tag(vtkfilename, nLocBas);
}

Interface_Element_Group::Interface_Element_Group(
  const std::string &vtkfilename, const int &nLocBas, const Vector_3 &centroid_in, const std::vector<double> &intervals_in) :
  groups_type {1}, axis_direction {-1}, centroid {centroid_in},
  intervals {intervals_in}, interval_tag {std::vector<int> {}}
{
  // Check the intervals
  SYS_T::print_fatal_if(VEC_T::get_size(intervals) < 2,
    "Error, Interface_Element_Group: illegal 'interval' vector. (1)\n");

  const int n_interval {VEC_T::get_size(intervals) - 1};
  
  for(int ii{0}; ii<n_interval; ++ii)
  {
    SYS_T::print_fatal_if(intervals[ii] >= intervals[ii + 1],
      "Error, Interface_Element_Group: illegal 'interval' vector. (2)\n");
  }

  SYS_T::print_fatal_if(intervals[0] != 0.0,
    "Error, Interface_Element_Group: illegal 'interval' vector. (3)");

  Tag(vtkfilename, nLocBas);
}

void Interface_Element_Group::Tag(const std::string &vtkfilename, const int &nLocBas)
{
  int s_nPts, s_nElem;
  std::vector<int> s_vecIEN;
  std::vector<double> s_ctrlPts;

  VTK_T::read_grid(vtkfilename, s_nPts, s_nElem, s_ctrlPts, s_vecIEN);

  interval_tag.resize(s_nElem);

  for(int ee{0}; ee<s_nElem; ++ee)
  {
    Vector_3 ele_centroid(0.0, 0.0, 0.0);
    for(int ii{0}; ii<nLocBas; ++ii)
    {
      int node {s_vecIEN[ee * nLocBas + ii]};
      ele_centroid.x() += s_ctrlPts[3 * node];
      ele_centroid.y() += s_ctrlPts[3 * node + 1];
      ele_centroid.z() += s_ctrlPts[3 * node + 2];
    }

    ele_centroid *= (1.0 / (double) nLocBas);

    interval_tag[ee] = Group(ele_centroid);
  }
}

int Interface_Element_Group::Group(const Vector_3 &ele_centroid)
{
  const int n_interval {VEC_T::get_size(intervals) - 1};

  switch (groups_type)
  {
    case 0:
    {
      for(int ii{0}; ii<n_interval; ++ii)
      {
        if(ele_centroid(axis_direction)>=intervals[ii] && ele_centroid(axis_direction)<=intervals[ii+1])
          return ii;

        SYS_T::print_fatal_if(ii==n_interval && ele_centroid(axis_direction)>intervals[ii+1],
          "Error, need better intervals end.\n");
      }
    }
    break;

    case 1:
    {
      const double dist {(ele_centroid - centroid).norm2()};
      for(int ii{0}; ii<n_interval; ++ii)
      {
        if(dist>=intervals[ii] && dist<=intervals[ii+1])
          return ii;

        SYS_T::print_fatal_if(ii==n_interval && dist>intervals[ii+1],
          "Error, need better intervals end.\n");
      }
    }
    break;
    
    default:
      SYS_T::print_fatal("Error, wrong group type.\n");
      break;
  }
}