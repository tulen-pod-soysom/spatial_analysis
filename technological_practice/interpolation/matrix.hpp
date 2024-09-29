#include <cassert>
#include <vector>

template <typename T> struct matrix {
  matrix() : _rows(0), _cols(0) {};
  matrix(int i, int j, T init = 0)
      : _rows(i), _cols(j), _data(std::vector<T>(i * j, init)) {};

  template <typename T_> matrix(std::vector<T_> v, int width) {
    assert(v % width != 0); // vector doesn't provide a rectangle matrix
    _data.resize(v.size());

    std::copy(v.begin(), v.end(), _data.begin());

    _rows = v.size() / width;
    _cols = width;
  }

  int rows() { return _rows; }
  int cols() { return _cols; }
  auto data() { return _data.data(); }

  
  struct column_iterator {
    int rows;
    int cols;
    int size;

    int col;
    int pos;

    void operator+=(int step) { pos += }
  };

protected:
  int _rows;
  int _cols;

  std::vector<T> _data;
};



template <typename Array, typename T>
class matrix_wrapper
{
  Array& arr;


  int _rows;
  int _cols;

  public:

  template <typename T_, size_t size_i, size_t size_j>
  matrix_wrapper<T_[size_i][size_j], T_>(T_ arr[size_i][size_j])
      : arr(arr), _rows(size_i), _cols(size_j){};

  matrix_wrapper<std::vector<std::vector<T>>, T>(
      std::vector<std::vector<T>> m)
      : arr(m), _rows(m.size()), _cols(m[0].size()){};

  

  auto rows() { return _rows; }
  auto cols() { return _cols; }

  T &operator()(int i, int j) { return arr[i][j]; }
};