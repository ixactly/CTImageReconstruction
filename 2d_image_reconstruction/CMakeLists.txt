cmake_minimum_required(VERSION 3.22)
project(2d_image_reconstruction)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenCV REQUIRED)

add_executable(2d_image_reconstruction main.cpp)
target_link_libraries(2d_image_reconstruction ${OpenCV_LIBS})
