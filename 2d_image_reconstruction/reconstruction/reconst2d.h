//
// Created by 森智希 on 2022/05/20.
//
#pragma once

#include <vector>
#include <array>
#include <cmath>

void ParallelBackProj(std::vector<float> &x_img, const std::vector<float> &b_proj);

void ParallelForwardProj(const std::vector<float> &x_img, std::vector<float> &b_proj);

void SIRT(std::vector<float> &x_img, const std::vector<float> &b_proj, const double alpha, const int num_iter);

void ART(std::vector<float> &x_img, const std::vector<float> &b_proj);

void Normalize(std::vector<float> &vec, const float max);
