cmake_minimum_required(VERSION 3.22)
add_library(ReconUtil INTERFACE MLEM.h Utils.h Params.h MLEM.h)
target_include_directories(ReconUtil INTERFACE ./)