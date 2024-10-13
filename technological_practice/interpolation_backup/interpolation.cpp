#include <iostream>
#include <boost/numeric/ublas/blas.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <armadillo>

template<typename T>
T linear_interpolation(T l , T r, T mu)
{
  return l + mu * (r - l);
}

template<typename X_type, typename Y_type>
Y_type linear_interpolation(X_type point, X_type lx , Y_type ly, X_type rx , Y_type ry)
{
  return linear_interpolation(ly,ry,(point - lx) / (rx-lx));
  // return ly + (ry-ly)/(rx-lx) * (point - lx);
}


template <typename Matrix, typename T = typename Matrix::value_type>
T det3(Matrix m)
{
  return 
  m(0,0) * (m(1,1)*m(2,2) - m(1,2)*m(2,1)) -
  m(1,0) * (m(1,0)*m(2,2) - m(1,2)*m(2,0)) +
  m(2,0) * (m(1,0)*m(1,1) - m(1,1)*m(2,0)) ;
}

template <typename Matrix, typename Vector>
Vector solve33_cramer(Matrix mat, Vector vec)
{
  auto d = det3(mat);

  auto mat1 = mat;
  auto mat2 = mat;
  auto mat3 = mat;

  for (auto i = 0; i <3; ++i) {
    mat1(i,0) = vec(i);
    mat2(i,1) = vec(i);
    mat3(i,2) = vec(i);
  }

  auto d1 = det3(mat1);
  auto d2 = det3(mat2);
  auto d3 = det3(mat3);
  // Vector out();
  // out(0) = d1/d;
  // out(1) = d2/d;
  // out(2) = d3/d;

  return {d1/d,d2/d,d3/d};
}

// template<typename X_type = double, typename Y_type = double>
// Y_type quadratic_interpolation(X_type x, Y_type p1_y, Y_type p2_y, Y_type p3_y, X_type p1_x, X_type p2_x,X_type p3_x)
// {
//   boost::numeric::ublas::matrix<X_type> mat(3, 3);
//   mat(0,0) = mat(1,0) = mat(2,0) = 1;
//   mat(0, 1) = p1_x;
//   mat(1, 1) = p2_x;
//   mat(2, 1) = p3_x;
//   mat(0, 2) = p1_x * p1_x;
//   mat(1, 2) = p2_x * p2_x;
//   mat(2, 2) = p3_x * p3_x;

//   // mat =
//   // {
//   //   1, p1_x, p1_x * p1_x,
//   //   1, p2_x, p2_x * p2_x,
//   //   1, p3_x, p3_x * p3_x,
//   // };

//   boost::numeric::ublas::vector<Y_type> v(3);
//   v(0) = p1_y;
//   v(1) = p2_y;
//   v(2) = p3_y;

//   auto out = solve33_cramer(mat, v);

//   return out(0) + out(1) * x + out(2) * x * x;
// }


template<typename X_type = double, typename Y_type = double>
Y_type quadratic_interpolation(X_type x, Y_type p1_y, Y_type p2_y, Y_type p3_y, X_type p1_x, X_type p2_x,X_type p3_x)
{
  arma::mat33 m = {
      {1, p1_x, p1_x * p1_x},
      {1, p2_x, p2_x * p2_x}, 
      {1, p3_x, p3_x * p3_x},
  };

  arma::colvec3 v = {p1_y,p2_y,p3_y};

  // auto out = solve33_cramer(m, v);
  arma::colvec3 out = arma::solve(m,v);
  return out(0) + out(1) * x + out(2) * x * x;
}

// template<typename X_type = double, typename Y_type = double>
// Y_type quadratic_interpolation(X_type x, Y_type p1_y, Y_type p2_y, Y_type p3_y, X_type p1_x, X_type p2_x,X_type p3_x)
// {
//   Y_type a[3];

//   a[2] = (p3_y - p1_y)/(p3_x - p1_x)/(p3_x-p2_x) - (p2_y - p1_y)/(p2_x - p1_x)/(p3_x-p2_x);
//   a[1] = (p2_y - p1_y)/(p2_x - p1_x) - a[2] *(p2_x + p1_x);
//   a[0] = p1_y - a[1]* p1_x - a[2]*p1_x*p1_x;


//   // std::cout << a[2] << ' ' << a[1] << ' ' << a[0] << std::endl;
//   return a[0] + a[1]*x + a[2]*x*x;
// }



//https://ru.wikipedia.org/wiki/Билинейная_интерполяция
template<typename ImageView, typename Image>
Image bilinear_interpolation(ImageView& m, size_t w, size_t h, size_t w_new, size_t h_new)
{
  if (w == w_new && h == h_new) return m;
  Image out(h_new,w_new,0);

  // у границ нет соседей
  for (auto i = 1; i < h_new - 1; ++i) {
    for (auto j = 1; j < w_new - 1; ++j)
    {
      double x = j/double(w_new);
      double y = i/double(h_new);

      int closest_lx = floor(x*(w));
      int closest_rx = ceil(x*(w));
      int closest_ly = floor(y*(h));
      int closest_ry = ceil(y*(h));

      if ((closest_lx == closest_rx) || (closest_ly == closest_ry)) {
        out(i, j) = m(closest_lx, closest_lx);
        continue;
      }

      auto p11 = m(closest_ly, closest_lx);
      auto p12 = m(closest_ly, closest_rx);
      auto p21 = m(closest_ry, closest_lx);
      auto p22 = m(closest_ry, closest_rx);

      if (closest_lx == closest_rx) {
        out(i, j) = (closest_ry - y * h) * p11 + (y * h - closest_ly) * p21;
        continue;
      }
      if (closest_ly == closest_ry) {
        out(i, j) = (closest_rx - x * w) * p11 + (x * w - closest_lx) * p12;
        continue;
      }

      auto fR1 = (closest_rx - x * w) * p11 + (x * w - closest_lx) * p12;
      auto fR2 = (closest_rx - x * w) * p21 + (x * w - closest_lx) * p22;

      out(i, j) = (closest_ry - y * h) * fR1 + (y * h - closest_ly) * fR2;
    }
  }

  return out;
}

// template<typename Image>
// Image bilinear_interpolation(Image& m, size_t w_new, size_t h_new)
// {
//     return bilinear_interpolation<typename Image::view_t,Image>(
//         boost::gil::view(m),m.width(),m.height(),w_new,h_new);
// }

int main()
{
  std::vector<double> x = {-1,0,1};
  std::vector<double> y = {1,0,1};

  for (auto a = -1.0; a <= 1.0; a+=0.2) {
    std::cout << quadratic_interpolation(a,y[0], y[1], y[2], x[0], x[1], x[2]) << std::endl;
  }
}