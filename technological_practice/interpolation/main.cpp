#include <armadillo>
#include <boost/gil/image_view.hpp>
#include <boost/gil/typedefs.hpp>
#include <boost/gil/extension/io/bmp.hpp>
#include <boost/gil.hpp>

#include <cstdlib>
#include <fstream>

#include "bilinear_interpolation.hpp"
#include "biquadratic_interpolation.hpp"

#include <random>

template <typename T>
auto sqr(T a)
{
    return a * a;
}

template <typename T>
T sqr(std::complex<T> a)
{
    return a.imag() * a.imag() + a.real() * a.real();
}

template <typename Matrix1, typename Matrix2>
auto MSE(Matrix1& m1, Matrix2& m2, unsigned rows, unsigned cols)
{
    double error = 0.0;
    for (auto i = 0; i < rows; ++i)
        for (auto j = 0; j < cols; ++j)
    {
        error += sqr(m1(i,j) - m2(i,j));
    }

    return error /(double)rows /(double)cols;
}

template <typename T1, typename T2>
auto MSE(arma::Mat<T1>& m1, arma::Mat<T2>& m2)
{
    assert(m1.n_cols == m2.n_cols);
    assert(m1.n_rows == m2.n_rows);
    return MSE(m1,m2,m1.n_rows,m1.n_cols);
}

int main()
{
    double k_max = 4;

    std::string filename; std::cin >> filename;

    arma::uchar_mat im = gray_image_from_file(filename);

    auto w = im.n_rows;
    auto h = im.n_cols;

    auto bilinear = im;
    auto biquadratic = im;

    unsigned n_experiments = 7;

    std::random_device rd;
    std::uniform_real_distribution<double> dist(0,k_max);

    for (auto i =0; i < n_experiments; ++i)
    {
        unsigned new_w = w * dist(rd);
        unsigned new_h = h * dist(rd);

        bilinear = bilinear_interpolation<arma::uchar_mat>(bilinear, bilinear.n_rows, bilinear.n_cols , new_w, new_h);
        biquadratic = biquadratic_interpolation<arma::uchar_mat>(biquadratic, biquadratic.n_rows, biquadratic.n_cols , new_w, new_h);
    }

    bilinear = bilinear_interpolation<arma::uchar_mat>(bilinear, bilinear.n_rows, bilinear.n_cols , w, h);
    biquadratic = biquadratic_interpolation<arma::uchar_mat>(biquadratic, biquadratic.n_rows, biquadratic.n_cols ,w,h);

    std::cout << "bilinear MSE: " << MSE(im,bilinear) << std::endl;
    std::cout << "biquadratic MSE: " << MSE(im,biquadratic) << std::endl;
}
