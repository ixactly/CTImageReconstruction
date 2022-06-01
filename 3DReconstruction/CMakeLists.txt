cmake_minimum_required(VERSION 3.22)
project(3DReconstruction)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenCV REQUIRED)

add_subdirectory(./reconstruction)

add_executable(3DReconstruction main.cpp)
target_link_libraries(3DReconstruction PRIVATE ${OpenCV_LIBS} ReconUtil)