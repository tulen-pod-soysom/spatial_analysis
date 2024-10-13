
#include <cstddef>
#include <iostream>
template <typename Array, typename T>
class matrix_wrapper
{
  Array& arr;


  int _rows;
  int _cols;

  public:

    template<size_t size_i, size_t size_j>
  matrix_wrapper<T[size_i][size_j], T>(T (&arr)[size_i][size_j])
      : arr(arr), _rows(size_i), _cols(size_j){};

//   matrix_wrapper<std::vector<std::vector<T>>, T>(
//       std::vector<std::vector<T>> m)
//       : arr(m), _rows(m.size()), _cols(m[0].size()){};

  

  auto rows() { return _rows; }
  auto cols() { return _cols; }

  T &operator()(int i, int j) { return arr[i][j]; }
};


template<typename Matrix>
auto print_matrix(Matrix m)
{
    for (auto i = 0; i != m.rows(); ++i) {
        for (auto j =0; j != m.cols(); ++j)
        {
            std::cout << m(i,j) << ' ';
        }
        std::cout << std::endl;
    }
}

int main()
{
    int arr[3][3] = {
        1,2,3,
        4,5,6,
        7,8,9,
    };

    matrix_wrapper<int[3][3],int> m(arr);

    print_matrix(m);

    m(2,2) += 3.1415926535;

    print_matrix(m);
}