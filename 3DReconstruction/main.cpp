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
    // ground truth
    Volume<float> ctGT(NUM_VOXEL, NUM_VOXEL, 1);
    Volume<float> ct(NUM_VOXEL, NUM_VOXEL, 1);

    Geometry geom(SRC_DETECT_DISTANCE, SRC_OBJ_DISTANCE, DETECTOR_SIZE);
    sinogram.load("../volume_bin/Reslice_of_yukiphantom_1024x1000.raw", NUM_DETECT_U, 1, NUM_PROJ);

    /*
    for (int i = NUM_VOXEL / 3; i < NUM_VOXEL * 2 / 3 + 1; i++) {
        for (int j = NUM_VOXEL / 3; j < NUM_VOXEL * 2 / 3 + 1; j++) {
            ctGT(i, j, 0) = 1.0;
        }
    }
     */

    MLEM<float> mlem;

    // measure clock
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    // main function
    // mlem.forwardproj(sinogram, ctGT, geom);
    mlem.reconstruct(sinogram, ct, geom, 1, 64);

    end = std::chrono::system_clock::now();
    double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / (1000.0 * 1000.0));
    std::cout << "time: " << time << " (s)" << std::endl;

    ct.show(0);

    /*
    std::string savefilePath =
            "../volume_bin/proj-" + std::to_string(NUM_DETECT_U) + "x" + std::to_string(NUM_DETECT_V) + "x" +
            std::to_string(NUM_PROJ) + ".raw";
    */
    std::string savefilePath =
            "../volume_bin/emos_yuki-" + std::to_string(NUM_VOXEL) + "x" + std::to_string(NUM_VOXEL) + ".raw";
    ct.save(savefilePath);
}