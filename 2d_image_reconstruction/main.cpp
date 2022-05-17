#include <fstream>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
#include <cmath>
#include <filesystem>
#include <algorithm>

inline const int H_IMG = 470;
inline const int W_IMG = 470;

inline const int NUM_DETECT = 470;
inline const int NUM_PROJ = 400;
inline const double PIXEL_SIZE = 0.8;
inline const double D_THETA = 2 * M_PI / NUM_PROJ;
inline const double AXIS_OFFSET = 12.3708265; // (unit: pixel) 12.3708265

void BackProjection(std::vector<float> &x_img, const std::vector<float> &b_proj);

void ForwardProjection(const std::vector<float> &x_img, std::vector<float> &b_proj);

void SIRT(std::vector<float> &x_img, const std::vector<float> &b_proj, const double alpha, const int num_iter);

void ART(std::vector<float> &x_img, const std::vector<float> &b_proj);

void Normalize(std::vector<float> &vec, const float max);

int main() {
    std::vector<float> x_image(H_IMG * W_IMG, 0);
    std::vector<float> b_proj(NUM_DETECT * NUM_PROJ);

    // load image binary
    std::filesystem::path cur_path = std::filesystem::current_path();
    std::ifstream ifile("../2d_images_bin/elephant-470x400.raw", std::ios::binary);
    if (!ifile) {
        std::cout << "file not opened" << std::endl;
        std::cout << cur_path.string() << std::endl;
        return 0;
    }
    ifile.read(reinterpret_cast<char *>(b_proj.data()), sizeof(float) * NUM_PROJ * NUM_DETECT);

    // --------------- main processing ---------------

    /* hexagonal
    for (int i = 0; i < H_IMG; ++i) {
        for (int j = 0; j < W_IMG; ++j) {
            if(std::fabs(i+j-H_IMG) < 100 && std::fabs(j - W_IMG/2) < 100) {
                x_image[H_IMG*i+j] = 1.5;
            }
        }
    }
    */
    // square

    for (int i = 0; i < H_IMG; ++i) {
        for (int j = 0; j < W_IMG; ++j) {
            if (std::fabs(i - H_IMG / 3) < 50 && std::fabs(j - W_IMG / 4) < 30) {
                x_image[H_IMG * i + j] = 2.0;
            }
        }
    }
    ForwardProjection(x_image, b_proj);
    std::fill(x_image.begin(), x_image.end(), 0);
    BackProjection(x_image, b_proj);
    // SIRT(x_image, b_proj, 0.001, 100);

    // --------------- end processing ---------------
    for (auto &e: x_image) {
        if (e < 0) {
            e = 0;
        }
    }

    Normalize(x_image, 1.0);

    cv::Mat img(x_image);
    cv::Mat prj(b_proj);

    // v = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]→
    // v.reshape(1, 3) = [1, 2, 3, 4,
    //                    5, 6, 7, 8,
    //                    9, 10, 11, 12]
    // 表示もこの行列の通りに表示される

    cv::Mat img_show = img.reshape(1, H_IMG);
    cv::Mat prj_show = prj.reshape(1, NUM_PROJ);
    cv::imshow("image", img_show);
    cv::imshow("sinogram", prj_show);

    cv::waitKey(0);

    std::ofstream ofs("../2d_images_bin/ele_image-470x470.raw", std::ios::binary);
    if (!ofs) {
        std::cout << "file not opened" << std::endl;
        std::cout << cur_path.string() << std::endl;
        return 0;
    }
    ofs.write(reinterpret_cast<char *>(x_image.data()), sizeof(float) * H_IMG * W_IMG);

    return 0;
}

void ART(std::vector<float> &x_img, const std::vector<float> &b_proj) {

    double theta = 0;
    const double x_rotate_center = (W_IMG - 1) / 2.0;
    const double y_rotate_center = (H_IMG - 1) / 2.0;
    double pos_ray_on_t;
    const double t_center = (NUM_DETECT - 1) / 2.0;

    // 本来はcudaなどの並列計算により, iterationを回すたびにシステム行列の要素を計算する．
    for (int k_proj = 0; k_proj < NUM_PROJ; ++k_proj) {
        // i, jは再構成画像のpixelのちょうど中心の座標を意味する.(not 格子点)

        // pixel driven forward projection
        std::vector<std::stack<std::pair<float, int>>> A_sparse(NUM_DETECT);
        for (int i_pic = 0; i_pic < H_IMG; ++i_pic) {
            for (int j_pic = 0; j_pic < W_IMG; ++j_pic) {
                pos_ray_on_t =
                        (j_pic - x_rotate_center) * std::cos(theta) + (i_pic - y_rotate_center) * std::sin(theta)
                        + t_center;
                const int t_floor = static_cast<int>(std::floor(pos_ray_on_t));

                if (pos_ray_on_t > 0.0 && pos_ray_on_t < NUM_DETECT - 1) {
                    // Linear interpolation
                    A_sparse[t_floor].push(std::make_pair(1 - (pos_ray_on_t - t_floor), W_IMG * i_pic + j_pic));
                    A_sparse[t_floor + 1].push(std::make_pair((pos_ray_on_t - t_floor), W_IMG * i_pic + j_pic));
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
    }
}

// projectionするセルの数で割る必要がありそう?
void BackProjection(std::vector<float> &x_img, const std::vector<float> &b_proj) {
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

void ForwardProjection(const std::vector<float> &x_img, std::vector<float> &b_proj) {
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
                    /*
                    A_sparse[t_floor].push(std::make_pair(1 - (pos_ray_on_t - t_floor), W_IMG * i_pic + j_pic));
                    A_sparse[t_floor + 1].push(std::make_pair((pos_ray_on_t - t_floor), W_IMG * i_pic + j_pic));
                    */
                } else if (pos_ray_on_t > -0.5 && pos_ray_on_t <= 0.0) { // corner case on using std::floor
                    b_proj[k_proj * NUM_DETECT + t_floor + 1] +=
                            (pos_ray_on_t - t_floor) * x_img[W_IMG * i_pic + j_pic];
                    // A_sparse[t_floor + 1].push(std::make_pair((pos_ray_on_t - t_floor), W_IMG * i_pic + j_pic));
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
        ForwardProjection(x_img, b_tmp);

        for (int i = 0; i < b_tmp.size(); ++i) {
            b_tmp[i] = alpha * (b_proj[i] - b_tmp[i]);
        }
        // fixing
        BackProjection(x_tmp, b_tmp);

        for (int i = 0; i < x_tmp.size(); ++i) {
            x_img[i] = x_img[i] + x_tmp[i];
        }
    }
}

void Normalize(std::vector<float> &vec, const float max) {
    float max_val = *std::max_element(vec.begin(), vec.end());
    for (auto &e: vec) e = e * max / max_val;
}