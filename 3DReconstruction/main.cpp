//
// Created by 森智希 on 2022/06/01.
//

#include <iostream>
#include <string>
#include "Volume.h"
#include "MLEM.h"
#include "Geometry.h"
#include "Params.h"
#include <chrono>
#include "Pbar.h"

int main() {

    Volume<float> sinogram(NUM_DETECT_U, NUM_DETECT_V, NUM_PROJ);
    // ground truth
    Volume<float> ct(NUM_VOXEL, NUM_VOXEL, NUM_VOXEL);
    Geometry geom(SRC_DETECT_DISTANCE, SRC_OBJ_DISTANCE, DETECTOR_SIZE);
    sinogram.load("../volume_bin/yukiphantom_float_1024x1024x1000.raw", NUM_DETECT_U, NUM_DETECT_V, NUM_PROJ);
    ct.load("../volume_bin/yuki_recon2-128x128x128.raw", NUM_VOXEL, NUM_VOXEL, NUM_VOXEL);
    /*
    for (int i = NUM_VOXEL / 3; i < NUM_VOXEL * 2 / 3 + 1; i++) {
        for (int j = NUM_VOXEL / 3; j < NUM_VOXEL * 2 / 3 + 1; j++) {
            for (int k = NUM_VOXEL / 3; k < NUM_VOXEL * 2 / 3 + 1; k++) {
                ct(i, j, k) = 1.0;
            }
        }
    }
    */

    // measure clock
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    // main function
    MLEM<float> mlem;
    // mlem.forwardproj(sinogram, ctGT, geom, Rotate::CCW);
    mlem.reconstruct(sinogram, ct, geom, 1, 50, Rotate::CCW, true);

    end = std::chrono::system_clock::now();
    double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() /
                                      (1000.0 * 1000.0));
    std::cout << "\n time: " << time << " (s)" << std::endl;

    // ct.show(NUM_VOXEL / 2);
    /*
    std::string savefilePath =
            "../volume_bin/cube_phantom-" + std::to_string(NUM_DETECT_U) + "x" + std::to_string(NUM_DETECT_V) + "x" +
            std::to_string(NUM_PROJ) + ".raw";
    */

    std::string savefilePath =
            "../volume_bin/yuki_recon3-" + std::to_string(NUM_VOXEL) + "x" +
            std::to_string(NUM_VOXEL) + "x" + std::to_string(NUM_VOXEL) + ".raw";

    ct.save(savefilePath);
}