clang++ main.cpp -o bilinear.out -O3 -larmadillo -DBILINEAR_INTERPOLATION
clang++ main.cpp -o biquadratic.out -O3 -larmadillo -DBIQUADRATIC_INTERPOLATION

./bilinear.out < cfg.txt
./biquadratic.out < cfg.txt

./plt
