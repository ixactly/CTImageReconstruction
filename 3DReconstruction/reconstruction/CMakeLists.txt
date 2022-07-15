cmake_minimum_required(VERSION 3.22)
add_library(reconst INTERFACE MLEM.h)

target_include_directories(reconst INTERFACE ../util)
