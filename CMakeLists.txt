cmake_minimum_required(VERSION 3.10.0)
project(polynome-comparaison VERSION 0.1.0 LANGUAGES C CXX)

find_package(Eigen3 REQUIRED)
find_package(libcmaes REQUIRED)

include_directories(include)

add_executable(cmaes-impl
    src/cmaes-impl.cpp
    src/regulator.cpp
    src/evaluation.cpp
)
# target_include_directories(cmaes-impl PUBLIC libcmaes::cmaes)
target_link_libraries(cmaes-impl libcmaes::cmaes Eigen3::Eigen)
