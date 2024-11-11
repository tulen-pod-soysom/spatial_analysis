
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

template<typename X_type = double, typename Y_type = double>
Y_type quadratic_interpolation(X_type x,X_type p1_x, X_type p2_x,X_type p3_x, Y_type p1_y, Y_type p2_y, Y_type p3_y)
{
  Y_type a[3];

  a[2] = (p3_y - p1_y)/(p3_x - p1_x)/(p3_x-p2_x) - (p2_y - p1_y)/(p2_x - p1_x)/(p3_x-p2_x);
  a[1] = (p2_y - p1_y)/(p2_x - p1_x) - a[2] *(p2_x + p1_x);
  a[0] = p1_y - a[1]* p1_x - a[2]*p1_x*p1_x;

  return a[0] + a[1]*x + a[2]*x*x;
}




template<typename Image, typename ImageView>
Image biquadratic_interpolation(ImageView& source, unsigned w, unsigned h, unsigned w_new, unsigned h_new)
{
  Image out(w_new,h_new);

  auto width_gcd = std::gcd(w, w_new);
  auto height_gcd = std::gcd(h, h_new);

  auto x_step = w / width_gcd;
  auto y_step = h / height_gcd;

  for (unsigned i = 0; i < w; i += x_step)
  for (unsigned j = 0; j < h; j += y_step)
  {
    // *--
    // |xx 
    // |xx 
    // 
    // where * is exact pixel value from source
    //       - is pixel interpolated on x axis
    //       | is pixel interpolated on y axis
    //       x is pixel interpolated in both axis using bi-linear/quadratic/cubic etc interpolation methods


    // fill the non interpolated point
    unsigned i_ = w_new * i / w;
    unsigned j_ = h_new * j / h;

    out(i_,j_) = source(i,j);


    // interpolate missing pixels
    // these ones are interpolated both in x and y axis
    for (auto i_interp = i_ + 1; (i_interp < i_ + w_new/width_gcd) && (i_interp < w_new); ++i_interp) {
      for (auto j_interp = j_ + 1;(j_interp < j_ + h_new/height_gcd) && (j_interp < h_new); ++j_interp)
      {
        double x = double(i_interp) * w / w_new;
        double y = double(j_interp) * h / h_new;

        double k1,k2,k3;
        double x1,x2,x3;
        double y1,y2,y3;
        x1 = floor(x) - 1;
        x2 = floor(x);
        x3 = floor(x) + 1;
        y1 = floor(y) - 1;
        y2 = floor(y);
        y3 = floor(y) + 1;

        if (floor(x) -1 == -1)
        {
            x1 = floor(x) + 0;
            x2 = floor(x) + 1;
            x3 = floor(x) + 2;
        }
        if (floor(y) -1 == -1)
        {
            y1 = floor(y) + 0;
            y2 = floor(y) + 1;
            y3 = floor(y) + 2;
        }

        if (floor(x) +1 == w)
        {
            x1 = floor(x) - 2;
            x2 = floor(x) - 1;
            x3 = floor(x) - 0;
        }
        if (floor(y) +1 == h)
        {
            y1 = floor(y) - 2;
            y2 = floor(y) - 1;
            y3 = floor(y) - 0;
        }



        // if (
        //     ((floor(x) + 1 >= w) || floor(x)-1 < 0) ||
        //     ((floor(y) + 1  >= h) || floor(y)-1 < 0)
        // )   continue;


        k1 = quadratic_interpolation<double,double>(
            x, 
            x1,x2,x3,
            source(x1,y1),
            source(x2,y1),
            source(x3,y1)
        );

        k2 = quadratic_interpolation<double,double>(
            x, 
            x1,x2,x3,
            source(x1,y2),
            source(x2,y2),
            source(x3,y2)
        );

        k3 = quadratic_interpolation<double,double>(
            x, 
            x1,x2,x3,
            source(x1,y3),
            source(x2,y3),
            source(x3,y3)
        );


        double z = quadratic_interpolation<double,double>(y,y1,y2,y3,k1,k2,k3);

        // converting from double to unsigned can be wrong
        out(i_interp, j_interp) = std::clamp(z,
            (double)std::numeric_limits<typeof(out(i, j))>::min(),
            (double)std::numeric_limits<typeof(out(i, j))>::max()
            );
      }
    }

    // these are interpolated only in x axis
    for (auto i_interp = i_ + 1; (i_interp < i_ + w_new/width_gcd) && (i_interp < w_new); ++i_interp) {
        double x = double(i_interp) * w / w_new;
        double x1,x2,x3;
        double z1,z2,z3;

        x1 = floor(x) - 1;
        x2 = floor(x);
        x3 = floor(x) + 1;

        if (floor(x) -1 == -1)
        {
            x1 = floor(x) + 0;
            x2 = floor(x) + 1;
            x3 = floor(x) + 2;
        }
        if (floor(x) +1 == w)
        {
            x1 = floor(x) - 2;
            x2 = floor(x) - 1;
            x3 = floor(x) - 0;
        }


        z1 = source(x1,j);
        z2 = source(x2,j);
        z3 = source(x3,j);

        double z = quadratic_interpolation<double,double>(x,x1,x2,x3,z1,z2,z3);

        // converting from double to unsigned can be wrong
        out(i_interp,j_) = std::clamp(z,
            (double)std::numeric_limits<typeof(out(i, j))>::min(),
            (double)std::numeric_limits<typeof(out(i, j))>::max()
            );
    }

    // these are interpolated only in y axis
    for (auto j_interp = j_ + 1; (j_interp < j_ + h_new/height_gcd) && (j_interp < h_new); ++j_interp) {
        double y = double(j_interp) * h / h_new;
        double y1,y2,y3;
        double z1,z2,z3;

        y1 = floor(y) - 1;
        y2 = floor(y);
        y3 = floor(y) + 1;

        if (floor(y) -1 == -1)
        {
            y1 = floor(y) + 0;
            y2 = floor(y) + 1;
            y3 = floor(y) + 2;
        }
        if (floor(y) +1 == h)
        {
            y1 = floor(y) - 2;
            y2 = floor(y) - 1;
            y3 = floor(y) - 0;
        }

        z1 = source(i,y1);
        z2 = source(i,y2);
        z3 = source(i,y3);


        double z = quadratic_interpolation<double,double>(y,y1,y2,y3,z1,z2,z3);

        // converting from double to unsigned can be wrong
        out(i_, j_interp) = std::clamp(z,
            (double)std::numeric_limits<typeof(out(i, j))>::min(),
            (double)std::numeric_limits<typeof(out(i, j))>::max()
            );
    }
  }

  return out;
}
