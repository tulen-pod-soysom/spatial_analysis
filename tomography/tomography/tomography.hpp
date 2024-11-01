#ifndef TOMOGRAPHY_HPP
#define TOMOGRAPHY_HPP

#include <cmath>
#include <cstdlib>
#include <utility>
#include <vector>

template<typename Func>
auto brezenhem_algorithm(int x0,int y0, int x1, int y1, Func f)
{

  int dx = std::abs(x1 - x0);
  int dy = std::abs(y1 - y0);
  int err = 0;
  int d_err;
  int dir;


  if (dx >= dy)
  {
    if (x1 < x0) {std::swap(x1,x0); std::swap(y1,y0);}
    d_err = (dy + 1);
    dir = (y1 - y0 > 0)? 1: -1;

    int y = y0;
    for (int i = 0; i < dx; ++i)
    {
        f(x0 + i, y);
        err += d_err;
        if (err >= dx +1)
        {
            y += dir;
            // y++;
            err -= dx + 1;
        }
    }
  }
  else
  {
    if (y1 < y0) {std::swap(y1,y0); std::swap(x1,x0);};
    d_err = (dx + 1);
    dir = (x1 - x0 > 0)? 1: -1;

    int x = x0;
    for (int j = 0; j < dy; ++j)
    {
        f(x,y0 + j);
        err += d_err;
        if (err >= dy + 1)
        {
            x += dir;
            // x++;
            err -= dy + 1;
        }
    }
  }
}

template<typename T>
auto in_range(T x, T x_min, T x_max)
{
    return (x >= x_min) && (x <= x_max);
}


template <typename Real = double>
auto line_rect_intersection(Real A, Real B, Real C, Real x0, Real y0, Real x1, Real y1)
{
    if (std::abs(A) < std::numeric_limits<Real>::epsilon()) {
        if (in_range(-C/B,y0,y1)) {return std::vector{x0,-C/B,x1,-C/B};}
    }

    if (std::abs(B) < std::numeric_limits<Real>::epsilon()) {
        if (in_range(-C/A,x0,x1)) {return std::vector{-C/A,y0,-C/A,y1};}
    }


    if (x1 < x0) std::swap(x1,x0);
    if (y1 < y0) std::swap(y1,y0);

    Real x_res_1;
    Real y_res_1;
    Real x_res_2;
    Real y_res_2;

    int found = 0;
    if (in_range(-(A*x0 + C) / B, y0, y1))
    {
        auto& x = x_res_1;
        auto& y = y_res_1;

        x = x0;
        y = -(A*x0 + C) / B;
        found++;
    }
    if (in_range(-(B*y0 + C)/ A, x0, x1))
    {
        auto& x = (found == 0)? x_res_1 : x_res_2;
        auto& y = (found == 0)? y_res_1 : y_res_2;

        x = -(B*y0 + C)/ A;
        y = y0;
        found++;
    }   
    
    if (in_range(-(A*x1 + C) / B, y0, y1))
    {
        auto& x = (found == 0)? x_res_1 : x_res_2;
        auto& y = (found == 0)? y_res_1 : y_res_2;
        x = x1;
        y = -(A*x1 + C) / B;
        found++;
    }
    if (in_range(-(B*y1 + C)/ A, x0, x1))
    {
        auto& x = (found == 0)? x_res_1 : x_res_2;
        auto& y = (found == 0)? y_res_1 : y_res_2;
        x = -(B*y1 + C)/ A;
        y = y1;
        found++;
    }

    if (found == 2)
    return std::vector{x_res_1,y_res_1,x_res_2,y_res_2};
    if (found == 1)
    return std::vector{x_res_1,y_res_1};
    else 
    return std::vector<Real>{};
}

template<typename Image>
auto get_projection(Image object, unsigned w, unsigned h, double angle = 0.0, double x_shift = 0)
{
    double A,B,C;
    A = cos(angle);
    B = sin(angle);

    C = -x_shift; // simplified according to trigonometric identity
    C += -A*double(w)/2.0 - B*double(h)/2.0;

    auto points = line_rect_intersection<double>(A, B, C, 0, 0, w-1, h-1);
    if (points.size() != 4) return 0.0;

    int p[4] = {
        int(floor(points[0])),
        int(floor(points[1])),
        int(floor(points[2])),
        int(floor(points[3]))
    };

    double sum = 0;
    brezenhem_algorithm(p[0], p[1], p[2], p[3], [&](int i, int j){sum += object(i,j);});

    return sum;
}


#endif //TOMOGRAPHY_HPP
