#include "image_io.hpp"
#include "tomography.hpp"
#include <armadillo>

int main()
{
    // auto im = gray_image_from_file<arma::mat>("image.png");

    // auto str = std::ofstream("test.txt");
    // im.print(str);

    // gray_image_to_file(im, im.n_rows, im.n_cols, "out.jpg");

    arma::umat m(30,30);
    brezenhem_algorithm(10, 10, 10, 10, [&](int x, int y) {m(x,y) = 1;});

    std::ofstream f("test.txt");
    m.print(f);
}
