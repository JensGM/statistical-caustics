cmake_minimum_required(VERSION 3.14)
project(ocean LANGUAGES C CXX)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)

find_package(FFTW3 REQUIRED)
find_package(GLM REQUIRED)

add_subdirectory(lib)
add_subdirectory(bin)
