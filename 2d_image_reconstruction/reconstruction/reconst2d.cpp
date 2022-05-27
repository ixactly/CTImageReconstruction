//
// Created by 森智希 on 2022/05/20.
//
#include <iostream>
#include "reconst2d.h"
#include "params.h"

// projectionするセルの数で割る必要がありそう?
void ParallelBackProj(const std::vector<float> &b_proj, std::vector<float> &x_img) {
    double theta = 0;
    const double x_rotate_center = (W_IMG - 1) / 2.0;
    const double y_rotate_center = (H_IMG - 1) / 2.0;
    double pos_ray_on_t;
    const double t_center = (NUM_DETECT - 1) / 2.0;

    // 本来はcudaなどの並列計算により, iterationを回すたびにシステム行列の要素を計算する．
    for (int k_proj = 0; k_proj < NUM_PROJ; ++k_proj) {
        // i, jは再構成画像のpixelのちょうど中心の座標を意味する.(not 格子点)
        // pixel driven back projection
        for (int i_pic = 0; i_pic < H_IMG; ++i_pic) {
            for (int j_pic = 0; j_pic < W_IMG; ++j_pic) {
                pos_ray_on_t =
                        (j_pic - x_rotate_center) * std::cos(theta) + (i_pic - y_rotate_center) * std::sin(theta)
                        + t_center;
                const int t_floor = static_cast<int>(std::floor(pos_ray_on_t));

                float adder_to_img = 0;
                if (pos_ray_on_t > 0.0 && pos_ray_on_t < NUM_DETECT - 1) {
                    // Linear interpolation
                    adder_to_img = static_cast<float>(
                            (b_proj[NUM_DETECT * k_proj + t_floor] * (t_floor + 1 - pos_ray_on_t) +
                             b_proj[NUM_DETECT * k_proj + t_floor + 1] * (pos_ray_on_t - t_floor)));
                } else if (pos_ray_on_t > -0.5 && pos_ray_on_t <= 0.0) { // corner case on using std::floor
                    adder_to_img = static_cast<float>(
                            (b_proj[NUM_DETECT * k_proj + t_floor + 1] * (pos_ray_on_t - t_floor)));
                } else if (pos_ray_on_t >= NUM_DETECT - 1 && pos_ray_on_t < NUM_DETECT - 1 + 0.5) {
                    adder_to_img = static_cast<float>(
                            (b_proj[NUM_DETECT * k_proj + t_floor] * (t_floor + 1 - pos_ray_on_t)));
                }

                x_img[W_IMG * i_pic + j_pic] += adder_to_img / NUM_PROJ;
            }
        }
        theta += D_THETA;
    }
}

void ParallelForwardProj(const std::vector<float> &x_img, std::vector<float> &b_proj) {
    // projection 0 fill
    std::fill(b_proj.begin(), b_proj.end(), 0);

    double theta = 0;
    const double x_rotate_center = (W_IMG - 1) / 2.0;
    const double y_rotate_center = (H_IMG - 1) / 2.0;
    double pos_ray_on_t;
    const double t_center = (NUM_DETECT - 1) / 2.0;

    // 本来はcudaなどの並列計算により, iterationを回すたびにシステム行列の要素を計算する．
    for (int k_proj = 0; k_proj < NUM_PROJ; ++k_proj) {
        // i, jは再構成画像のpixelのちょうど中心の座標を意味する.(not 格子点)
        // pixel driven back projection
        // std::vector<std::stack<std::pair<float, int>>> A_sparse(NUM_DETECT);

        for (int i_pic = 0; i_pic < H_IMG; ++i_pic) {
            for (int j_pic = 0; j_pic < W_IMG; ++j_pic) {
                pos_ray_on_t =
                        (j_pic - x_rotate_center) * std::cos(theta) + (i_pic - y_rotate_center) * std::sin(theta)
                        + t_center;
                const int t_floor = static_cast<int>(std::floor(pos_ray_on_t));

                if (pos_ray_on_t > 0.0 && pos_ray_on_t < NUM_DETECT - 1) {
                    // Linear interpolation
                    b_proj[k_proj * NUM_DETECT + t_floor] +=
                            (1 - (pos_ray_on_t - t_floor)) * x_img[W_IMG * i_pic + j_pic];
                    b_proj[k_proj * NUM_DETECT + t_floor + 1] +=
                            (pos_ray_on_t - t_floor) * x_img[W_IMG * i_pic + j_pic];

                } else if (pos_ray_on_t > -0.5 && pos_ray_on_t <= 0.0) { // corner case on using std::floor
                    b_proj[k_proj * NUM_DETECT + t_floor + 1] +=
                            (pos_ray_on_t - t_floor) * x_img[W_IMG * i_pic + j_pic];

                } else if (pos_ray_on_t >= NUM_DETECT - 1 && pos_ray_on_t < NUM_DETECT - 1 + 0.5) {
                    b_proj[k_proj * NUM_DETECT + t_floor] +=
                            (1 - (pos_ray_on_t - t_floor)) * x_img[W_IMG * i_pic + j_pic];
                }
            }
        }
        theta += D_THETA;
    }
}

void SIRT(std::vector<float> &x_img, const std::vector<float> &b_proj, const double alpha, const int num_iter) {
    std::vector<float> b_tmp(b_proj.size());
    std::vector<float> x_tmp(x_img.size());
    const float max_b = *std::max_element(b_proj.begin(), b_proj.end());

    for (int iter = 0; iter < num_iter; ++iter) {
        // calculate error
        ParallelForwardProj(x_img, b_tmp);

        for (int i = 0; i < b_tmp.size(); ++i) {
            b_tmp[i] = alpha * (b_proj[i] - b_tmp[i]);
        }
        // fixing
        ParallelBackProj(b_tmp, x_tmp);

        for (int i = 0; i < x_tmp.size(); ++i) {
            x_img[i] = x_img[i] + x_tmp[i];
        }
    }
}

void Fan2Para(const std::vector<float> &fan_proj, std::vector<float> &par_proj) { // Geometry classの追加
    const double L0 = DISTANCE_SOURCE_SAMPLE;
    const double Ld = DISTANCE_SOURCE_SAMPLE + DISTANCE_DETECT_SAMPLE;

    double t_fan, theta_fan;
    for (int t_par = 0; t_par < NUM_DETECT; ++t_par) {
        double Xp = (t_par - NUM_DETECT / 2.0) * PIXEL_SIZE;
        double Xp_max = NUM_DETECT * PIXEL_SIZE / 2.0;

        double scale = std::sqrt(L0 * L0 - Xp_max * Xp_max) / Ld;
        t_fan = scale * Ld * Xp / std::sqrt(L0 * L0 - Xp * Xp) / PIXEL_SIZE + NUM_DETECT / 2.0;
        if (t_fan < 0 || t_fan > NUM_DETECT) {
            std::cout << "t: " << t_fan << " theta: " << theta_fan << std::endl;
            continue;
        }

        for (int theta_par = 0; theta_par < NUM_PROJ; ++theta_par) {
            theta_fan =
                    std::fmod(2 * M_PI + theta_par * D_THETA - std::atan2(Xp, std::sqrt(L0 * L0 - Xp * Xp)), 2 * M_PI) /
                    D_THETA;
            if (theta_fan < 0 || theta_fan > NUM_PROJ) {
                std::cout << "t: " << t_fan << " theta: " << theta_fan << std::endl;
                continue;
            }

            // lerp(Linear Interpolation) in t, theta direction
            double a = t_fan - std::floor(t_fan);
            double b = theta_fan - std::floor(theta_fan);
            int floor_idx = (static_cast<int>(theta_fan) + 1) % NUM_PROJ;
            double w1 = fan_proj[std::floor(theta_fan) * NUM_DETECT + std::floor(t_fan)] * (1 - a) * (1 - b);
            double w2 = fan_proj[std::floor(theta_fan) * NUM_DETECT + (std::floor(t_fan) + 1)] * a * (1 - b);
            double w3 = fan_proj[floor_idx * NUM_DETECT + std::floor(t_fan) + 1] * a * b;
            double w4 = fan_proj[floor_idx * NUM_DETECT + std::floor(t_fan)] * (1 - a) * b;

            par_proj[theta_par * NUM_DETECT + t_par] = w1 + w2 + w3 + w4;
        }
    }
}

void ART(std::vector<float> &x_img, const std::vector<float> &b_proj) {
    /*double theta = 0;
    const double x_rotate_center = (W_IMG - 1) / 2.0;
    const double y_rotate_center = (H_IMG - 1) / 2.0;
    double pos_ray_on_t;
    const double t_center = (NUM_DETECT - 1) / 2.0;

    // 本来はcudaなどの並列計算により, iterationを回すたびにシステム行列の要素を計算する．
    for (int k_proj = 0; k_proj < NUM_PROJ; ++k_proj) {
        // i, jは再構成画像のpixelのちょうど中心の座標を意味する.(not 格子点)
        // Asparse
        // pixel driven forward projection
        for (int i_pic = 0; i_pic < H_IMG; ++i_pic) {
            for (int j_pic = 0; j_pic < W_IMG; ++j_pic) {

                pos_ray_on_t =
                        (j_pic - x_rotate_center) * std::cos(theta) + (i_pic - y_rotate_center) * std::sin(theta)
                        + t_center;
                const int t_floor = static_cast<int>(std::floor(pos_ray_on_t));

                float norm_a = 0;
                float dot_ax = 0;

                if (pos_ray_on_t > 0.0 && pos_ray_on_t < NUM_DETECT - 1) {
                    // Linear interpolation
                    norm_a += std::pow(1 - (pos_ray_on_t - t_floor), 2) +
                              std::pow((pos_ray_on_t - t_floor), 2);
                    dot_ax += (1 - (pos_ray_on_t - t_floor)) * x_img[W_IMG * i_pic + j_pic] +
                              (pos_ray_on_t - t_floor) * x_img[W_IMG * i_pic + j_pic];

//                    A_sparse[t_floor].push(std::make_pair(1 - (pos_ray_on_t - t_floor), W_IMG * i_pic + j_pic));
//                    A_sparse[t_floor + 1].push(std::make_pair((pos_ray_on_t - t_floor), W_IMG * i_pic + j_pic));

                } else if (pos_ray_on_t > -0.5 && pos_ray_on_t <= 0.0) { // corner case on using std::floor
                    A_sparse[t_floor + 1].push(std::make_pair((pos_ray_on_t - t_floor), W_IMG * i_pic + j_pic));
                } else if (pos_ray_on_t >= NUM_DETECT - 1 && pos_ray_on_t < NUM_DETECT - 1 + 0.5) {
                    A_sparse[t_floor].push(std::make_pair(1 - (pos_ray_on_t - t_floor), W_IMG * i_pic + j_pic));
                }
            }
        }

        // pop from stack data and calculate projection
        int t = 0;
        for (auto &e: A_sparse) {
            float norm_a = 0;
            float dot_ax = 0;

            std::vector<std::pair<float, int>> vec_a(e.size());
            for (auto &a: vec_a) {
                a = e.top();
                norm_a += std::pow(a.first, 2);
                dot_ax += a.first * x_img[a.second];
                e.pop();
            }

            for (auto &a: vec_a) {
                x_img[a.second] += (b_proj[k_proj * NUM_DETECT + t] - dot_ax) * a.first / norm_a;
            }
            t++;
        }
        theta += D_THETA;
    }*/
}

void Normalize(std::vector<float> &vec, const float max) {
    float max_val = *std::max_element(vec.begin(), vec.end());
    for (auto &e: vec) e = e * max / max_val;
}
