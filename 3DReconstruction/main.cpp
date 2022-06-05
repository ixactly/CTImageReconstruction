//
// Created by 森智希 on 2022/06/01.
//

#include <iostream>
#include <Volume.h>
#include <MLEM.h>
#include <Geometry.h>
#include <Params.h>

int main() {
    Volume<float> sinogram(NUM_DETECT_U, NUM_DETECT_V, NUM_PROJ);
    Volume<float> ct(NUM_VOXEL, NUM_VOXEL, NUM_VOXEL);
    Geometry geom(SRC_DETECT_DISTANCE, SRC_OBJ_DISTANCE, DETECTOR_SIZE);
    // sinogram.load("../volume_bin/sphere-tori-float-500x500x500.raw", 500, 500, 500);

    MLEM<float> mlem;
    mlem.reconstruct(sinogram, ct, geom);
    ct.show(NUM_VOXEL/2);




}