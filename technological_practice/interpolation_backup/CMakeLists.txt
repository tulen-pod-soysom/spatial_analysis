project(interpolation)

find_package(Armadillo REQUIRED)
link_libraries(armadillo)

add_executable(interpolation interpolation.cpp)