#include <armadillo>
#include <cstdlib>
#include <fstream>

#include "../bilinear_interpolation.hpp"
#include "../biquadratic_interpolation.hpp"

int main()
{
    arma::uchar_mat im(10,10);

    for (auto i =0; i < im.n_rows; ++i) {
        for (auto j =0 ; j < im.n_cols; ++j)
        {
            // im(i,j) = double(rand() )* 255/ RAND_MAX;
            im(i,j) = 255.0/2.0*(sin(2*3.1415926535*(i+j)/10.0)+1);
        }
    }

    unsigned width, height;
    std::cin >> width >> height;

#ifdef BILINEAR_INTERPOLATION
    arma::uchar_mat res = bilinear_interpolation<arma::uchar_mat>(im, im.n_cols, im.n_rows, width,height);
    std::ofstream f2("res_bilinear.txt");
#endif

#ifdef BIQUADRATIC_INTERPOLATION
    arma::uchar_mat res = biquadratic_interpolation<arma::uchar_mat>(im, im.n_cols, im.n_rows, width,height);
    std::ofstream f2("res_biquadratic.txt");
#endif

    std::ofstream f1("im.txt");

    im.print(f1);
    res.print(f2);
}
