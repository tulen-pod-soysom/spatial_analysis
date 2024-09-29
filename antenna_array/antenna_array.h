

#ifndef ANTENNA_ARRAY_H
#define ANTENNA_ARRAY_H

#include <cmath>
#include <cstddef>
#include <vector>
#include <tuple>
#include <complex>
#include <numeric>

const static double PI = 3.14159265358979323;

struct vec2{
    vec2(double x = 0, double y = 0): x(x), y(y){}
    double x{0}; double y{0};
};

struct antenna
{
    antenna(vec2 p= {0,0}): position(p){}
    antenna(double x, double y): position(x,y){}
    vec2 position{0,0};
    // std::complex<double> modifier{1};
};

int index_2d_to_1d(int i,int j, int width);
std::tuple<int,int> index_1d_to_2d(int a, int width);

double distance2(double x1, double y1, double x2, double y2);

double distance(double x1, double y1, double x2, double y2);

std::tuple<double,double> decart_to_polar(double x, double y);


// double polar_distance2(double r1, double theta1, double r2, double theta2)
// {

// }


struct antenna_array
{
    std::vector<antenna> antennas;
    std::vector<std::complex<double>> modifiers;

    double wave_length = 1.0;

    void create_default_configuration();

    void create_grid_configuration_from(std::vector<std::vector<bool>> a);

    template <typename InputIt>
    std::complex<double> get_output(InputIt begin, InputIt end)
    {
        return std::inner_product(begin,end,modifiers.begin(),std::complex<double>(0.0));
    }
};



// Описывает источник ЭМ сигнала
struct emitter
{
    //phi = 0 для нормали к решетке
    double R, phi, theta;

    double wave_length = 1;
    double frequency;
    double A = 1;


    // theta in degrees
    double radial_distance_to(double r, double phi_, double theta_ = 0)
    {
        const double degree_to_radians = PI/180.0;
        return sqrt(R*R + r*r - 2*R*r*cos(degree_to_radians*(phi - theta)));
    }

    double decart_distance_to(double x, double y, double z = 0)
    {
        const double degree_to_radians = PI/180.0;
        double x_ = R * cos(degree_to_radians*phi) * sin(degree_to_radians*theta);
        double y_ = R * sin(degree_to_radians*phi) * sin(degree_to_radians*theta);
        double z_ = R * cos(degree_to_radians*theta);


        return sqrt((x_-x)*(x_-x) + (y_-y)*(y_-y) + (z_-z)*(z_-z));
    }

    // Высчитывает значение комплексной амплитуды в точке
    std::complex<double> amplitude_at_point_radial(double r, double theta)
    {
        double r_ = radial_distance_to(r,theta);
        std::complex<double> ikr_(0.0,2*PI/wave_length * r_);
        return A*exp(-ikr_) / r_;
    }

    std::complex<double> amplitude_at_point_decart(double x, double y, double z = 0)
    {
        double r_ = decart_distance_to(x,y,z);
        std::complex<double> ikr_(0.0,2*PI/wave_length * r_);
        return A*exp(-ikr_) / r_;
    }

    std::complex<double> phase_at_point(double r, double theta)
    {
        double r_ = radial_distance_to(r, theta);
        return exp(std::complex<double>(0.0, 2 * PI / wave_length * r_));
    }
};

std::vector<double> get_directional_diagram(antenna_array& arr, size_t samples, double distance);

#endif
