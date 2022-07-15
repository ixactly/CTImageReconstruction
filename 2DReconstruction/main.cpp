#include <fstream>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <vector>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include <array>

#include "reconst2d.h"
#include "../3DReconstruction/util/Params.h"

int main() {
    std::vector<float> x_image(H_IMG * W_IMG, 0);
    std::vector<float> b_proj(NUM_DETECT * NUM_PROJ, 0);
    std::vector<float> parallel(NUM_DETECT * NUM_PROJ, 0);

    // load image binary
    std::filesystem::path cur_path = std::filesystem::current_path();
    std::ifstream ifile("../2d_images_bin/sphere-tori-float-500x500.raw", std::ios::binary);
    if (!ifile) {
        std::cout << "file not opened" << std::endl;
        std::cout << cur_path.string() << std::endl;
        return 0;
    }
    ifile.read(reinterpret_cast<char *>(b_proj.data()), sizeof(float) * NUM_PROJ * NUM_DETECT);

    // --------------- main processing ---------------
    /*
    for (int i = 0; i < H_IMG; ++i) {
        for (int j = 0; j < W_IMG; ++j) {
            if(std::fabs(i+j-H_IMG) < 100 && std::fabs(j - W_IMG/2) < 100) {
                x_image[H_IMG*i+j] = 1.0;
            }
        }
    }
     */


    // initialize
    for (int i = 0; i < H_IMG; ++i) {
        for (int j = 0; j < W_IMG; ++j) {
            if (std::sqrt(std::pow(i - H_IMG / 2, 2) + std::pow(j - W_IMG / 2, 2)) < H_IMG / 2)
            x_image[H_IMG * i + j] = 1.0;
        }
    }

    Fan2Para(b_proj, parallel);

    for (auto &e : parallel) {
        if (e < 0.0) {
            e = 1e-10;
        }
    }

    for (int iter = 0; iter < 50; iter++) MLEM(x_image, parallel);
    ParallelForwardProj(x_image, parallel);
    // --------------- end processing ---------------

    // Normalize(x_image, 1.0);

    cv::Mat img(x_image);
    cv::Mat prj(parallel);

    // v = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]→
    // v.reshape(1, 3) = [1, 2, 3, 4,
    //                    5, 6, 7, 8,
    //                    9, 10, 11, 12]
    // 表示もこの行列の通りに表示される

    cv::Mat img_show = img.reshape(1, H_IMG);
    cv::Mat prj_show = prj.reshape(1, NUM_PROJ);
    cv::imshow("img", img_show);
    cv::imshow("parallel", prj_show);

    cv::waitKey(0);

    std::ofstream ofs("../2d_images_bin/sphere-tori-image-mlem50-500x500.raw", std::ios::binary);
    if (!ofs) {
        std::cout << "file not opened" << std::endl;
        std::cout << cur_path.string() << std::endl;
        return 0;
    }

    ofs.write(reinterpret_cast<char *>(x_image.data()), sizeof(float) * H_IMG * H_IMG);

    return 0;
}

