#include <algorithm>
#include <boost/gil/extension/dynamic_image/any_image.hpp>
#include <boost/gil/extension/io/png/tags.hpp>
#include <boost/gil/image.hpp>
#include <boost/gil/image_view.hpp>
#include <boost/gil/io/base.hpp>
#include <boost/gil/io/conversion_policies.hpp>
#include <boost/gil/io/write_view.hpp>
#include <boost/gil/io/read_image_info.hpp>
#include <boost/gil/io/read_view.hpp>
#include <boost/gil/io/typedefs.hpp>
#include <boost/gil/typedefs.hpp>
#include <functional>
#include <numeric>
#include <random>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/gil.hpp>
#include <boost/gil/io/read_and_convert_image.hpp>
#include <boost/gil/io/read_image.hpp>
#include <boost/gil/extension/io/png.hpp>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>

using gray_image = boost::numeric::ublas::matrix<unsigned char>;
using image_spectre = boost::numeric::ublas::matrix<std::complex<double>>;

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

template <typename MatrixIn = gray_image, typename MatrixOut = image_spectre>
MatrixOut fft_2d(MatrixIn &m, bool inverse = false) {

  MatrixOut output(m.size1(), m.size2());

  auto row_out = output.begin1();
  for (auto row = m.begin1(); row != m.end1(); ++row, ++row_out) {
    fft(row.begin(), row.end(), row_out.begin(), inverse ? 1 : -1);
  }

  auto col_out = output.begin2();
  for (auto col = m.begin2(); col != m.end2(); ++col, ++col_out) {
    fft(col_out.begin(), col_out.end(), col_out.begin(), inverse ? 1 : -1);
  }

  return output;
}

template<typename T>
std::ostream &operator<<(std::ostream &str, boost::numeric::ublas::matrix<T> &im) {
  for (auto i = 0; i < im.size1(); i++) {
    for (auto j = 0; j < im.size2(); j++) {
      str << (int)im.at_element(i, j) << ' ';
    }
    str << std::endl;
  }
  return str;
}

template <typename Matrix>
void print_abs(std::ostream &str, Matrix &im) {
  for (auto i = 0; i < im.size1(); i++) {
    for (auto j = 0; j < im.size2(); j++) {
      str << abs(im.at_element(i, j)) << ' ';
    }
    str << std::endl;
  }
};

template <typename Matrix>
void print_log(std::ostream &str, Matrix &im) {
  for (auto i = 0; i < im.size1(); i++) {
    for (auto j = 0; j < im.size2(); j++) {
      str << log(1 + abs(im.at_element(i, j))) << ' ';
    }
    str << std::endl;
  }
};

template <typename Matrix>
void save_to_draw_spectre(std::ostream& str, Matrix& m)
{
  // constant part in unnecessary
  auto a = m(0,0);
  m(0,0) = 0;
  print_log(str, m);
  m(0,0) = a;
}


inline double sqr(double x) { return x * x; };
inline double sqr(std::complex<double> x) { return x.imag() * x.imag() + x.real() * x.real(); };

template <typename Matrix>
double calculate_energy(const Matrix& m)
{
  return std::accumulate(m.data().begin(),m.data().end(),0.0,[](auto a, auto b){return a + sqr(b);});
}


template <typename Matrix>
auto filter_image(Matrix& m, double energy_fraction = 0.95)
{
  auto total_energy = calculate_energy(m);
  total_energy -= sqr(m(0,0));

  int i = m.size1() / 2 - 1;
  int j = m.size2() / 2 - 1;

  double deceased_energy = 0;
  double threshold = total_energy * (1 - energy_fraction); 

  while (deceased_energy < threshold) {

    if (j > 0)
      for (auto i_ = 0; i_ < m.size1(); i_++) {
        deceased_energy += sqr(m(i_, j));
        deceased_energy += sqr(m(i_, m.size2() - 1 - j));

        if (deceased_energy > threshold)
          return;

        m(i_, j) = 0;
        m(i_, m.size2() - 1 - j) = 0;
      }

    if (i > 0)
      for (auto j_ = 0; j_ < m.size2(); j_++) {
        deceased_energy += sqr(m(i, j_));
        deceased_energy += sqr(m(m.size1() - 1 - i, j_));

        if (deceased_energy > threshold)
          return;

        m(i, j_) = 0;
        m(m.size1() - 1 - i, j_) = 0;
      }

    if (i > 0) i--;
    if (j > 0) j--;

    if (i <= 0 && j <= 0) return;
  }

};



struct gauss {
  double A, s_x, s_y, m_x, m_y;
  double operator()(double x, double y) {
    return A * exp(-(sqr((x - m_x) / s_y) + sqr((y - m_y) / s_y)));
  }
};

auto load_from_file(std::string filepath)
{
  using namespace boost::gil;

  rgb8_image_t image;
  read_image(filepath,image,png_tag());

  gray_image mat(image.height(),image.width());

  auto view = image._view;
  for (auto i = 0; i < mat.size1(); ++i) {
    for (auto j = 0; j < mat.size2(); ++j)
    {
      rgb8_pixel_t pixel = *view.at(j,i);
      double R = at_c<0>(pixel);
      double G = at_c<1>(pixel);
      double B = at_c<2>(pixel);
      mat(i,j) = 0.299* R + 0.587 * G + 0.114 * B;
    }
  }
  return mat;
}

template<typename Matrix, typename A = typename Matrix::value_type,typename T = int, typename Preprocessor = std::function<T(A)>>
auto save_to_file(std::string filename, Matrix& m, Preprocessor preprocessor = [](A a)-> int{return a;})
{
  boost::gil::rgb8_image_t image(m.size2(),m.size1());

  for (auto i = 0; i < m.size2(); ++i) {
    for (auto j = 0; j < m.size1(); ++j)
    {
      int v = std::clamp(preprocessor(m(j,i)),0,255);
      image._view(i,j).at(std::integral_constant<int, 0>()) = v;
      image._view(i,j).at(std::integral_constant<int, 1>()) = v;
      image._view(i,j).at(std::integral_constant<int, 2>()) = v;
    }
  }

  boost::gil::write_view(filename,image._view,boost::gil::png_tag());
}

template<typename Matrix>
auto save_spectre(std::string filename, Matrix& m, bool log = true)
{
  auto a = m(0,0);
  m(0,0) = 0;

  auto pair = std::minmax_element(m.data().begin(),m.data().end(),[](auto a, auto b){return sqr(a) < sqr(b);});

  double max = abs(*pair.second);
  double max2 = sqr(*pair.second);
  if (log)
    save_to_file(filename, m, [&max2](auto a)-> int{return 255*log2(1 + sqr(a)) /log2(1 + max2);});
  else
    save_to_file(filename, m, [&max](auto a)-> int{return 255*abs(a)/max;});
  m(0,0) = a;
}


auto generate_image()
{
  gauss g1 = {63, 128, 128, 127, 127};
  gauss g2 = {63, 64, 64, 63, 256};
  gauss g3 = {63, 64, 64, 256, 256};

  gray_image im(512, 512, 127);
  for (auto i = 0; i < im.size1(); i++) {
    for (auto j = 0; j < im.size2(); j++) {
      im.at_element(i, j) += g1(i, j) + g2(i, j) + g3(i, j);
    }
  }
  return im; 
}

auto noise_image(size_t width, size_t height)
{
  static std::random_device rd;
  std::normal_distribution<double> dist;

  boost::numeric::ublas::matrix<double> im(width,height);

  for (auto& i : im.data()) {
      i = dist(rd);
  }

  return im;
}

//returns actual noise energy fraction
template<typename Matrix>
double add_noise(Matrix& m, double noise_energy_fraction = 0.05)
{
  auto noise = noise_image(m.size1(), m.size2());

  double image_energy = calculate_energy(m);

  double noise_energy = calculate_energy(noise);

  noise *= sqrt(noise_energy_fraction * image_energy / noise_energy);

  for (auto i = 0; i < m.size1(); ++i) {
    for (auto j = 0; j < m.size2(); ++j)
    {
      int a = int(m(i,j)) + noise(i,j);
      m(i,j) = std::clamp(a,0,255);
    }
  }

  return calculate_energy(m) / image_energy - 1;
}


//https://ru.wikipedia.org/wiki/Билинейная_интерполяция
template<typename Matrix>
Matrix bilinear_interpolation(const Matrix& m, size_t w, size_t h, size_t w_new, size_t h_new)
{
  if (w == w_new && h == h_new) return m;
  Matrix out(h_new,w_new,0);

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

template <typename T>
boost::numeric::ublas::matrix<T> bilinear_interpolation(boost::numeric::ublas::matrix<T>& m, unsigned w_new, unsigned h_new)
{
  return bilinear_interpolation(m,m.size2(),m.size1(),w_new,h_new);
}


unsigned closest_pow_2(int a)
{
  if ((a & (a - 1)) == 0)
    return a;
  return pow(2,ceil(log2(a)));
}


template<typename Matrix>
Matrix interpolate_to_pow_2(const Matrix& m)
{
  auto max_ = std::max(closest_pow_2(m.size2()),closest_pow_2(m.size1()));
  return bilinear_interpolation(m, m.size2(), m.size1(), max_, max_);
}

template <typename T>
void ask(std::string question, T& answer)
{
  std::cout << question << std::endl;
  std::cin >> answer;
}

template <>
void ask(std::string question, bool& answer)
{
  std::cout << question << std::endl;
  char ans;
  for (;;) {
    std::cin >> ans;
    if (ans == 'y' || ans == 'Y') {
      answer = true;
      return;
    } else if (ans == 'n' || ans == 'N') {
      answer = false;
      return;
    }
  }
}

int main() {

  bool log;
  bool from_file;
  std::string filepath;
  double noise_energy_fraction;
  double leave_energy_fraction;

  loading:

  ask("from file?", from_file);
  if (from_file) ask("filepath: ", filepath);
  gray_image im = from_file? load_from_file(filepath) : generate_image();
  im = interpolate_to_pow_2(im);

  adding_noise:
  ask("noise energy fraction: ", noise_energy_fraction);
  auto noised_im = im;
  auto noise_fraction = add_noise(noised_im,noise_energy_fraction);
  // std::cout << "\tactual noise energy fraction = " << noise_fraction << std::endl;
  // std::cout << "\tto filter - " << noise_fraction / (1 + noise_fraction) << std::endl;

  save_to_file("output/image1.png", noised_im);


  ask("log scale? ", log);

  image_spectre spec = fft_2d(noised_im);
  save_spectre("output/image2.png", spec);
  
  filtering:
  ask("leave energy fraction: ", leave_energy_fraction);
  image_spectre filtered_spec = fft_2d(noised_im);
  filter_image(filtered_spec,leave_energy_fraction);
  save_spectre("output/image3.png", filtered_spec);

  image_spectre im_restored = fft_2d(filtered_spec, true);
  save_to_file("output/image4.png", im_restored,[](auto a){return int(abs(a));});

  int answer;
  ask("What to do next:\n 1 - generate new image,\n 2 - set new noise,\n 3 - set new filter threshold,\n 4 - quit", answer);

  switch (answer) {
    case 1: goto loading;
    case 2: goto adding_noise;
    case 3: goto filtering;
    case 4: return 0;
  }
  
}
