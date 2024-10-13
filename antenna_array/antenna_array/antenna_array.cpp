#include <antenna_array.h>
#include <tuple>

int index_2d_to_1d(int i, int j, int width)
{
    return i*width + j;
}

std::tuple<int,int> index_1d_to_2d(int a, int width)
{
    return std::make_tuple(a / width,a % width);
}

double distance2(double x1, double y1, double x2, double y2)
{
    double dx = x2-x1;
    double dy = y2-y1;
    return dx*dx + dy*dy;
}

double distance(double x1, double y1, double x2, double y2)
{return sqrt(distance2(x1,y1,x2,y2));}


auto from_sphere_coords(double r, double phi, double theta = 0)
{
    double x = r * cos(phi) * sin(theta);
    double y = r * sin(phi) * sin(theta);
    double z = r * cos(theta);
    return std::make_tuple(x,y,z);
}


std::tuple<double, double> decart_to_polar(double x, double y)
{
    return std::make_tuple(sqrt(x*x+y*y),180*atan2(y,x)/PI);
}

std::vector<double> get_directional_diagram(antenna_array &arr, size_t width,
                                            double distance) {
  emitter em;
  em.R = distance;

  std::vector<double> diagram(width * width);

  for (auto i = 0; i < width; ++i) {
    em.phi = -180 + i / (double)width * 360;

    for (auto j = 0; j < width; ++j) {
      em.theta = 90 - j / (double)width * 180;

      std::vector<std::complex<double>> amplitudes(arr.antennas.size());

      for (auto k = 0; k < arr.antennas.size(); ++k) {
        auto [x, y] = arr.antennas[k].position;
        // auto [r,phi] = decart_to_polar(x,y);

        // amplitudes[k] = em.amplitude_at_point_radial(r,phi);
        amplitudes[k] = em.amplitude_at_point_decart(x, y);
      }
      diagram[i * width + j] =
          abs(arr.get_output(amplitudes.begin(), amplitudes.end()));
    }
  }
  return diagram;
}

std::vector<double> get_directional_diagram_projection(antenna_array &arr, size_t width,
                                            double distance) {
    emitter em;
    em.R = distance;



    std::vector<double> diagram(width * width);


    for (auto i = 0; i < width; ++i) {

        double x = -distance + i * 2*distance / (double) width;

        // em.phi = -180 + i / (double)width * 360;

        for (auto j = 0; j < width; ++j) {
            // em.theta = 90 - j / (double)width * 180;
            double y = -distance + j * 2*distance / (double) width;
            double z = sqrt(distance*distance - x*x - y*y);

            if (x*x+y*y > distance*distance) continue;
            em.phi = atan2(y,x) * 180. / PI;
            em.theta = acos(z/distance) * 180. / PI;

            std::vector<std::complex<double>> amplitudes(arr.antennas.size());

            for (auto k = 0; k < arr.antennas.size(); ++k) {
                auto [x_, y_] = arr.antennas[k].position;
                // auto [r,phi] = decart_to_polar(x,y);

                // amplitudes[k] = em.amplitude_at_point_radial(r,phi);
                amplitudes[k] = em.amplitude_at_point_decart(x_, y_);
            }
            diagram[i * width + j] =
                abs(arr.get_output(amplitudes.begin(), amplitudes.end()));
        }
    }
    return diagram;
}


void antenna_array::create_default_configuration(double wave_factor_period)
{
    size_t size_x = 10, size_y = 10;
    double w = size_x * wave_length*wave_factor_period, h = size_y * wave_length*wave_factor_period;

    antennas = std::vector<antenna>(size_x*size_y);
    modifiers = std::vector<std::complex<double>>(size_x*size_y,1.0);

    for (auto a = 0; a < size_x*size_y; ++a)
    {
        auto [i,j] = index_1d_to_2d(a,size_x);
        antennas[a].position = vec2(-w/2.0 + j*w/size_x,-h/2.0 + i*h/size_y);
    }
}

void antenna_array::create_grid_configuration_from(std::vector<std::vector<bool> > a,double wave_factor_period)
{
    size_t size_x = a[0].size(), size_y = a.size();
    double w = size_x * wave_length*wave_factor_period, h = size_y * wave_length*wave_factor_period;

    if (a.empty()) return;

    size_t size = 0;
    for (int i = 0; i < a.size(); ++i)
        for (int j = 0; j < a[0].size(); ++j)
            size += a[i][j]? 1:0;

    antennas.clear();
    antennas.reserve(size);
    modifiers = std::vector<std::complex<double>>(size,1.0);


    for (int i = 0; i < a.size(); ++i)
        for (int j = 0; j < a[0].size(); ++j)
            if (a[i][j])
            {
                antennas.emplace_back(antenna(-w/2.0 + j*w/size_x,-h/2.0 + i*h/size_y));
            }
}
