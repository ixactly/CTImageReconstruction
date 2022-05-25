#include <fstream>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <array>

#include "reconst2d.h"
#include "params.h"

int main() {
    std::vector<float> x_image(H_IMG * W_IMG, 0);
    std::vector<float> b_proj(NUM_DETECT * NUM_PROJ);
    std::vector<float> parallel(NUM_DETECT * NUM_PROJ);

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
    // ParallelForwardProj(x_image, b_proj);
    std::fill(x_image.begin(), x_image.end(), 0);
    Fan2Para(b_proj, parallel);
    ParallelBackProj(parallel, x_image);

    // SIRT(x_image, b_proj, 0.001, 100);

    // --------------- end processing ---------------
    for (auto &e: x_image) {
        if (e < 0) {
            e = 0;
        }
    }

    Normalize(x_image, 1.0);

    cv::Mat img(b_proj);
    cv::Mat prj(parallel);

    // v = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]→
    // v.reshape(1, 3) = [1, 2, 3, 4,
    //                    5, 6, 7, 8,
    //                    9, 10, 11, 12]
    // 表示もこの行列の通りに表示される

    cv::Mat img_show = img.reshape(1, NUM_PROJ);
    cv::Mat prj_show = prj.reshape(1, NUM_PROJ);
    cv::imshow("img", img_show);
    cv::imshow("parallel", prj_show);

    cv::waitKey(0);

    std::ofstream ofs("../2d_images_bin/para-470x400.raw", std::ios::binary);
    if (!ofs) {
        std::cout << "file not opened" << std::endl;
        std::cout << cur_path.string() << std::endl;
        return 0;
    }
    // ofs.write(reinterpret_cast<char *>(parallel.data()), sizeof(float) * NUM_PROJ * NUM_DETECT);

    return 0;
}

