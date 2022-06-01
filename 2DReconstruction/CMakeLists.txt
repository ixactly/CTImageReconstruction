cmake_minimum_required(VERSION 3.22)
project(2DReconstruction)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenCV REQUIRED)

add_subdirectory(./reconstruction)

add_executable(2DReconstruction main.cpp)
target_link_libraries(2DReconstruction PRIVATE ${OpenCV_LIBS} reconst2d)
