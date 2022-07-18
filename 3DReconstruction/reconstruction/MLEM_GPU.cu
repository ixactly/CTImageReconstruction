//
// Created by tomokimori on 22/07/16.
//

#include "MLEM_GPU.cuh"
#include "Geometry_GPU.cuh"
#include "../util/Volume.h"

void MLEM_CUDA::reconstruct(Volume<float> &sinogram, Volume<float> &voxel, const Geometry &geom, const int epoch, const int batch, Rotate dir) {

}

__global__ void forwardProj(const Geometry& geom) {
    const unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;

    auto [u, v] = geom.vox2det(x, y, z, rot * n, sinogram.size(), voxel.size());
    if (geom.isHitDetect(u, v, sinogram.size())) {
        lerpScatter(u, v, n, projTmp, voxel(x, y, z));
    }
}

__global__ void backwardProj() {

}