cmake_minimum_required(VERSION 3.22)
project(3DReconstruction)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenCV REQUIRED)

# openmp
if (APPLE)
    set(CMAKE_CXX_FLAGS "-fopenmp -L/opt/homebrew/Cellar/llvm/13.0.1_1/lib -Wl,-rpath,/opt/homebrew/Cellar/llvm/13.0.1_1/lib")
endif (APPLE)

if (WIN32)
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

if (UNIX AND NOT APPLE)
endif ()

add_subdirectory(./reconstruction)

add_executable(3DReconstruction main.cpp)
target_link_libraries(3DReconstruction PRIVATE ${OpenCV_LIBS} reconst)