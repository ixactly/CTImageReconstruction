cmake_minimum_required(VERSION 3.22)
add_library(ReconUtil INTERFACE MLEM.cpp Utils.h Params.h MLEM.h)
target_include_directories(ReconUtil INTERFACE ./)