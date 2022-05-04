#include <fstream>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
#include <cmath>
#include <filesystem>

inline const int H_IMG = 470;
inline const int W_IMG = 470;

inline const int NUM_DETECT = 470;
inline const int NUM_PROJ = 400;
inline const double PIXEL_SIZE = 0.8;
inline const double D_THETA = 2 * M_PI / NUM_PROJ;
inline const double AXIS_OFFSET = 12.3708265; // (unit: pixel)

void BackProjection(std::vector<float> &x_img, const std::vector<float> &b_proj);

int main() {
    std::vector<float> x_image(H_IMG * W_IMG);
    std::vector<float> b_proj(NUM_DETECT * NUM_PROJ);

    // load image binary
    std::filesystem::path cur_path = std::filesystem::current_path();
    std::ifstream file("../2d_images_bin/elephant-470x400.raw", std::ios::binary);
    if (!file) {
        std::cout << "file not opened" << std::endl;
        std::cout << cur_path.string() << std::endl;
    }
    file.read(reinterpret_cast<char *>(b_proj.data()), sizeof(float) * NUM_PROJ * NUM_DETECT);


    BackProjection(x_image, b_proj);

    cv::Mat img(x_image);

    // v = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]→
    // v.reshape(1, 3) = [1, 2, 3, 4,
    //                    5, 6, 7, 8,
    //                    9, 10, 11, 12]

    cv::Mat img_show = img.reshape(1, H_IMG);

    cv::imshow("reconstruction", img_show);
    cv::waitKey(0);
    return 0;
}

void Art2D(std::vector<std::vector<float>> &A, const std::vector<float> &b) {
    // check if matrix can be multiplied
    /*
    if (!(A.size() == b.size() && A[0].size() == x.size())) {
        std::cout << "invalid matrix size" << std::endl;
        return;
    }
     */
}

void ForwardProjection(std::vector<float> &x_img, std::vector<float> &b_proj) {

}

void BackProjection(std::vector<float> &x_img, const std::vector<float> &b_proj) {
    double theta = 0;
    const double x_rotate_center = (W_IMG - 1) / 2.0 + AXIS_OFFSET;
    const double y_rotate_center = (H_IMG - 1) / 2.0;
    double pos_ray_on_t;
    const double t_center = (NUM_DETECT - 1) / 2.0;

    // debug variable
    const int check_k = 100;
    // 本来はcudaなどの並列計算により, iterationを回すたびにシステム行列の要素を計算する．
    for (int k_proj = 0; k_proj < NUM_PROJ; ++k_proj) {
        // i, jは再構成画像のpixelのちょうど中心の座標を意味する.(not 格子点)
        // pixel driven back projection
        for (int i_pic = 0; i_pic < H_IMG; ++i_pic) {
            for (int j_pic = 0; j_pic < W_IMG; ++j_pic) {
                pos_ray_on_t =
                        (j_pic - x_rotate_center) * std::cos(theta) + (i_pic - y_rotate_center) * std::sin(theta) +
                        t_center;
                const int t_floor = static_cast<int>(std::floor(pos_ray_on_t));

                if (k_proj == check_k)
                    std::cout << pos_ray_on_t << " ";

                if (pos_ray_on_t > 0.0 && pos_ray_on_t < NUM_DETECT - 1) {
                    // Linear interpolation
                    x_img[W_IMG * i_pic + j_pic] = b_proj[NUM_PROJ * i_pic + t_floor] * (t_floor + 1 - pos_ray_on_t) +
                                                   b_proj[NUM_PROJ * i_pic + t_floor + 1] * (pos_ray_on_t - t_floor);
                } else if (pos_ray_on_t > -0.5 && pos_ray_on_t <= 0.0) { // corner case on using std::floor
                    x_img[W_IMG * i_pic + j_pic] = b_proj[NUM_PROJ * i_pic + t_floor + 1] * (pos_ray_on_t - t_floor);
                } else if (pos_ray_on_t >= NUM_DETECT - 1 && pos_ray_on_t < NUM_DETECT - 1 + 0.5) {
                    x_img[W_IMG * i_pic + j_pic] = b_proj[NUM_PROJ * i_pic + t_floor] * (t_floor + 1 - pos_ray_on_t);
                }
            }
            if (k_proj == check_k)
                std::cout << std::endl;
        }
        theta += D_THETA;
    }
}

void FilteredBackProjection2D() {

}