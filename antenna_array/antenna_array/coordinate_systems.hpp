// decart are assumed as default coordinate system

#include <cmath>
#include <tuple>
constexpr double PI = 3.1415265358979323;

template <typename Real = double>
auto from_sphere_coords(Real r, Real phi, Real theta = PI/2.0)
{
    Real x = r * cos(phi) * sin(theta);
    Real y = r * sin(phi) * sin(theta);
    Real z = r * cos(theta);
    return std::make_tuple(x,y,z);
}

template <typename Real = double>
auto to_sphere_coords(Real x, Real y, Real z)
{
    Real r = sqrt(x*x + y*y + z*z);
    Real phi = atan2(y, x); 
    Real theta = std::acos(z/r);

    return std::make_tuple(r,phi,theta);
}

template <typename Real = double>
auto from_cylindric_coords(Real r, Real phi, Real z)
{
    Real x = r * cos(phi);
    Real y = r * sin(phi);
    return std::make_tuple(x,y,z);
}

template <typename Real = double>
auto to_cylindric_coords(Real x, Real y, Real z)
{
    Real r = sqrt(x*x + y*y);
    Real phi = atan2(y, x); 

    return std::make_tuple(r,phi,z);
}


template <typename Real = double>
auto from_bipolar_coords(Real r, Real alpha, Real beta)
{
    Real x = 
    return std::make_tuple(x,y,z);
}

template <typename Real = double>
auto to_bipolar_coords(Real x, Real y, Real z)
{
    Real r = sqrt(x*x + y*y + z*z);
    Real alpha = atan2(z, x); 
    Real beta = atan2(z, y); 

    return std::make_tuple(r,alpha,beta);
}


template <typename Real = double>
auto rotate_OX(Real angle, Real x, Real y, Real z)
{
    Real x_ = x;
    Real y_ = y*cos(angle) - z*sin(angle);
    Real z_ = y*sin(angle) + z*cos(angle);
    return std::make_tuple(x_,y_,z_);
}


template <typename Real = double>
auto rotate_OY(Real angle, Real x, Real y, Real z)
{
    Real x_ = x*cos(angle) + z*sin(angle);
    Real y_ = y;
    Real z_ = -x*sin(angle) + z*cos(angle);
    return std::make_tuple(x_,y_,z_);
}

template <typename Real = double>
auto rotate_OZ(Real angle, Real x, Real y, Real z)
{
    Real x_ = x*cos(angle) - y*sin(angle);
    Real y_ = -x*sin(angle) + y*cos(angle);
    Real z_ = z;
    return std::make_tuple(x_,y_,z_);
}

