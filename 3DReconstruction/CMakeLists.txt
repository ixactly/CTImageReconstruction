cmake_minimum_required(VERSION 3.22)
project(3DReconstruction)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenCV REQUIRED)

# openmp
if (APPLE)
    set(CMAKE_CXX_FLAGS "-fopenmp -L/opt/homebrew/Cellar/llvm/13.0.1_1/lib -Wl,-rpath,/opt/homebrew/Cellar/llvm/13.0.1_1/lib")
endif (APPLE)

if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

if (UNIX AND NOT APPLE)
    enable_language(CUDA)
    find_package(OpenMP REQUIRED)
    find_package(CUDA REQUIRED)
    if(CUDA_FOUND)
        message("Cuda version: " ${CUDA_VERSION})
    endif()
endif ()

add_subdirectory(reconstruction)
add_subdirectory(util)

add_executable(3DRecon)
target_sources(3DRecon
        PRIVATE
        main.cpp)

target_include_directories(3DRecon PRIVATE reconstruction util)
target_link_libraries(3DRecon PRIVATE ${OpenCV_LIBS} reconst util)