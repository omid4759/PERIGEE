#ifndef QUADPTS_USERDEFINED_TRIANGLE_HPP
#define QUADPTS_USERDEFINED_TRIANGLE_HPP

#include "Vec_Tools.hpp"
#include "IQuadPts.hpp"

class QuadPts_UserDefined_Triangle : public IQuadPts
{
  public:
    QuadPts_UserDefined_Triangle( const int &in_num_pts = 1 ) : num_pts( in_num_pts )
    {
      qp = std::vector<double> (3 * num_pts, 0.333333333333333);
      qw = std::vector<double> (num_pts, 0.5 / num_pts);
    }
    
    virtual ~QuadPts_UserDefined_Triangle() = default;

    virtual int get_dim() const {return 3;}

    virtual int get_num_quadPts() const {return num_pts;}

    virtual double get_qp(unsigned int ii, unsigned int comp) const
    {return qp[3*ii+comp];}

    virtual double get_qw(unsigned int ii) const
    {return qw[ii];}

    virtual void set_qp(unsigned int ii, unsigned int comp, const double &value) const
    {qp[3*ii+comp] = value;}

    virtual void set_qw(unsigned int ii, const double &value) const
    {qw[ii] = value;}

  private:
    const int num_pts;

    // qp : length 3 x num_pts. Stores the r-s-t coordinates of the
    //      quadrature points.
    //      t = 1 - r - s.
    // qw : length num_pts. Stores the quadrature weights.
    std::vector<double> qp {};
    std::vector<double> qw {};
};

#endif