//
// Created by 森智希 on 2022/06/01.
//

#include <iostream>
#include <string>
#include <Volume.h>
#include <MLEM.h>
#include <Geometry.h>
#include <Params.h>
#include <chrono>

int main() {

    Volume<float> sinogram(NUM_DETECT_U, NUM_DETECT_V, NUM_PROJ);
    Volume<float> ct(NUM_VOXEL, NUM_VOXEL, NUM_VOXEL);
    Geometry geom(SRC_DETECT_DISTANCE, SRC_OBJ_DISTANCE, DETECTOR_SIZE);
    // sinogram.load("../volume_bin/sphere-tori-float-500x500x500.raw", 500, 500, 500);

    for (int i = NUM_VOXEL / 3; i < NUM_VOXEL * 2 / 3 + 1; i++) {
        for (int j = NUM_VOXEL / 3; j < NUM_VOXEL * 2 / 3 + 1; j++) {
            for (int k = NUM_VOXEL / 3; k < NUM_VOXEL * 2 / 3 + 1; k++) {
                ct(i, j, k) = 1.0;
            }
        }
    }

    MLEM<float> mlem;

    // measure clock
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    // main function
    mlem.forwardproj(sinogram, ct, geom);


    end = std::chrono::system_clock::now();
    double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / (1000.0 * 1000.0));
    std::cout << time << std::endl;

    ct.show(NUM_VOXEL/2);
    // sinogram.show(0);

    std::string savefilePath =
            "../volume_bin/proj-" + std::to_string(NUM_DETECT_U) + "x" + std::to_string(NUM_DETECT_V) + "x" +
            std::to_string(NUM_PROJ) + ".raw";
    sinogram.save(savefilePath);
}