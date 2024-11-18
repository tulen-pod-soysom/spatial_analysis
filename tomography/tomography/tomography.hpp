#ifndef TOMOGRAPHY_HPP
#define TOMOGRAPHY_HPP

#include <cmath>
#include <cstdlib>
#include <utility>
#include <vector>
#include <algorithm>
#include <armadillo>
#include <complex>
#include <bilinear_interpolation.hpp>


// just not to cause an ODR
template <typename T = int>
unsigned closest_pow_2(int a)
{
    if ((a & (a - 1)) == 0)
        return a;
    return pow(2,ceil(log2(a)));
}

template <typename InputIt, typename OutputIt>
void fft(InputIt begin, InputIt end, OutputIt out, double is) {
    int size = end - begin;
    if (size & size - 1)
        throw "size is not a power of two";
    static const double PI = 3.14159265358979323;

    // if (begin != out) {
    std::copy(begin, end, out);
    // }

    std::complex<double> temp, w, c;
    long i, j, istep;
    long m, mmax;
    long n = end - begin;
    double fn, r1, theta;
    fn = (double)n;
    double r = PI * is;

    j = 1;
    for (i = 1; i <= n; i++) {
        long i1 = i - 1;
        if (i < j) {
            int j1 = j - 1;
            std::swap(*(out + j1), *(out + i1));
        }
        m = n / 2;
        while (j > m) {
            j -= m;
            m = (m + 1) / 2;
        }
        j += m;
    }
    mmax = 1;
    while (mmax < n) {
        istep = 2 * mmax;
        r1 = r / (double)mmax;
        for (m = 1; m <= mmax; m++) {
            theta = r1 * (double)(m - 1);
            w = std::polar(1.0, theta);
            for (i = m - 1; i < n; i += istep) {
                j = i + mmax;
                c = *(out + j);
                temp = w * c;
                *(out + j) = *(out + i) - temp;
                *(out + i) = *(out + i) + temp;
            }
        }
        mmax = istep;
    }
    if (is > 0)
        for (i = 0; i < n; i++) {
            *(out + i) /= fn;
        }
};



template <typename T>
struct row_it : public std::iterator<std::random_access_iterator_tag,T>
{
    T& m;
    int i;
    int j;

    row_it(T& m, int i, int j):m(m),i(i),j(j){}
    auto operator+(int a)
    {
        return row_it<T>{m,i+a,j};
    }
    auto operator-(int a)
    {return operator+(-a);}

    auto& operator++()
    {i++; return *this;}

    auto& operator* ()
    {return m(i,j);}

    auto operator- (row_it<T> a){return i - a.i;}
};

template <typename T>
struct col_it: public std::iterator<std::random_access_iterator_tag,T>
{
    T& m;
    int i;
    int j;

    col_it(T& m, int i, int j):m(m),i(i),j(j){}
    auto operator+(int a)
    {
        return col_it<T>{m,i,j+a};
    }

    auto operator-(int a)
    {return operator+(-a);}

    auto& operator++()
    {j++; return *this;}

    auto& operator* ()
    {return m(i,j);}

    auto operator- (col_it a){return j - a.j;}
};

template <typename MatrixIn , typename MatrixOut>
MatrixOut fft_2d(MatrixIn &m, int rows, int cols, bool inverse = false) {

    MatrixOut output(rows, cols);

    for (auto i = 0; i < rows; ++i)
    {
        col_it<MatrixOut> col_out{output,i,0};
        col_it<MatrixIn> col    {m,i,0};
        col_it<MatrixIn> col_end{m,i,cols};

        fft(col,col_end,col_out,inverse ? 1:-1);
    }

    for (auto j = 0; j < cols; ++j)
    {
        row_it<MatrixOut> row_out{output,0,j};
        row_it<MatrixOut> row    {output,0,j};
        row_it<MatrixOut> row_end{output,rows,j};

        fft(row,row_end,row_out,inverse? 1:-1);
    }

    // auto row_out = output.begin1();
    // for (auto row = m.begin1(); row != m.end1(); ++row, ++row_out) {
        // fft(row.begin(), row.end(), row_out.begin(), inverse ? 1 : -1);
    // }

    // auto col_out = output.begin2();
    // for (auto col = m.begin2(); col != m.end2(); ++col, ++col_out) {
        // fft(col_out.begin(), col_out.end(), col_out.begin(), inverse ? 1 : -1);
    // }




    return output;
}



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

    if ((found == 2) || (found == 3) || (found == 4))
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
    if (points.size() != 4) 
        return 0.0;

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

template<typename ImageView, typename Image = ImageView>
auto sinogram_to_spectre(ImageView& sinogram, int w, int h, Image& output)
{
    std::complex<double> center_sum = 0;

    // sinogram = bilinear_interpolation<arma::mat>(sinogram,w,h,w,closest_pow_2(h));
    // h = closest_pow_2(h);

    // w means amount of divisions from 0 to pi
    for (auto i = 0; i < w; ++i)
    {
        double angle = i*3.14159265358979323 / (double)w;
        // center_sum+= sinogram(i,h/2+1);

        arma::cx_vec v(h);
        for (auto j = 0; j < h; ++j)
            v(j) = sinogram(i,j);

        fft(v.begin(),v.end(),v.begin(),-1);
        center_sum += v(0);

        for (auto j = 0; j < h/2; ++j)
            std::swap(v(j),v(j+h/2));

        for (auto j = 0; j < h; ++j)
        {
            double x = h/2.0 + cos(angle)*(j-h/2);
            double y = h/2.0 + sin(angle)*(j-h/2);

            x = round(x);
            y = round(y);

            x = std::clamp(x,0.0,double(h-1));
            y = std::clamp(y,0.0,double(h-1));

            // auto tmp = sinogram(i,j);
            // output(x,y) = tmp;
            output(x,y) = v(j);
        }
    }

    output(h/2+1,h/2+1) = center_sum / double (w);
}

template<typename ImageView, typename Image = ImageView>
auto restore_from_spectre(ImageView spectre, int w, int h, Image& output)
{

}


#endif //TOMOGRAPHY_HPP
