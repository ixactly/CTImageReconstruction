cmake_minimum_required(VERSION 3.22)

if (USE_CUDA_)
    set(CMAKE_CUDA_STANDARD 17)

    add_library(reconst SHARED)
    target_include_directories(reconst INTERFACE ../util)
    target_sources(reconst PRIVATE
            MLEM_GPU.cu
            MLEM_GPU.cuh
            Geometry_GPU.cuh)

else()

add_library(reconst INTERFACE MLEM.h)
target_include_directories(reconst INTERFACE ../util)

endif()
