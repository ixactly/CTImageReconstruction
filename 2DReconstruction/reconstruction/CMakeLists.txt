cmake_minimum_required(VERSION 3.22)
add_library(reconst2d STATIC reconst2d.cpp reconst2d.h params.h)
target_include_directories(reconst2d INTERFACE .)