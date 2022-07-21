//
// Created by tomokimori on 22/07/16.
//

#ifndef INC_3DRECONSTRUCTION_MLEM_GPU_CUH
#define INC_3DRECONSTRUCTION_MLEM_GPU_CUH

#include "../util/Volume.h"
#include "Geometry_GPU.cuh"

class MLEM_CUDA {
public :
    MLEM_CUDA() = default;

    void reconstruct(Volume<float> &sinogram, Volume<float> &voxel, const GeometryCUDA* geom, const int epoch, const int batch, Rotate dir);
};

#endif //INC_3DRECONSTRUCTION_MLEM_GPU_CUH
