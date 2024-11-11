
#include <algorithm>
#include <cmath>
#include <numeric>
template<typename T>
T linear_interpolation(T l , T r, T mu)
{
  return l + mu * (r - l);
}

template<typename X_type, typename Y_type>
Y_type linear_interpolation(X_type point, X_type lx, X_type rx , Y_type ly , Y_type ry)
{
  return linear_interpolation(ly,ry,(point - lx) / (rx-lx));
  // return ly + (ry-ly)/(rx-lx) * (point - lx);
}


template<typename Image, typename ImageView>
Image bilinear_interpolation(ImageView& source, unsigned w, unsigned h, unsigned w_new, unsigned h_new)
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

        double x1 = floor(x),x2 = ceil(x);
        double y1 = floor(y),y2 = ceil(y);

        if (ceil(x) >= w)
        {
          x1--;
          x2--;
        }
        if (ceil(y) >= h)
        {
          y1--;
          y2--;
        }

        double k1 = linear_interpolation<double,double>(
            x, x1,x2,
            source(x1,y1),
            source(x2,y1)
        );

        double k2 = linear_interpolation<double,double>(
            x, x1,x2,
            source(x1,y2),
            source(x2,y2)
        );

        double z = linear_interpolation<double,double>(y,y1,y2,k1,k2);
        out(i_interp,j_interp) = std::clamp(z,
            (double)std::numeric_limits<typeof(out(i, j))>::min(),
            (double)std::numeric_limits<typeof(out(i, j))>::max()
            );
      }
    }

    // these are interpolated only in x axis
    for (auto i_interp = i_ + 1; (i_interp < i_ + w_new/width_gcd) && (i_interp < w_new); ++i_interp) {
        double x = double(i_interp) * w / w_new;

        double x1 = floor(x),x2 = ceil(x);
        if (ceil(x) >= w)
        {
          x1--;
          x2--;
        }

        double z = linear_interpolation<double,double>(x,x1,x2,source(x1,j),source(x2,j));
        out(i_interp,j_) = std::clamp(z,
            (double)std::numeric_limits<typeof(out(i, j))>::min(),
            (double)std::numeric_limits<typeof(out(i, j))>::max()
            );
    }

    // these are interpolated only in y axis
    for (auto j_interp = j_ + 1; (j_interp < j_ + h_new/height_gcd) && (j_interp < h_new); ++j_interp) {
        double y = double(j_interp) * h / h_new;
        double y1 = floor(y),y2 = ceil(y);
        if (ceil(y) >= h)
        {
          y1--;
          y2--;
        }
        double z = linear_interpolation<double,double>(y,y1,y2,source(i,y1),source(i,y2));
        out(i_, j_interp) = std::clamp(z,
            (double)std::numeric_limits<typeof(out(i, j))>::min(),
            (double)std::numeric_limits<typeof(out(i, j))>::max()
            );
    }
  }

  return out;
}
