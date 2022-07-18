//
// Created by tomokimori on 22/07/16.
//

#include "MLEM_GPU.cuh"
#include "cuda.h"
#include "cuda_runtime.h"
#include "Geometry_GPU.cuh"

__global__ void xzPlaneFor(const int sizeD[3], const int sizeV[3], float* devSino, const float* devVoxel, GeometryCUDA& geom, const int y, const int n) {
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int z = blockIdx.y * blockDim.y + threadIdx.y;

    const int coord[4] = {x, y, z, n};

    geom.forwardProj(coord, sizeD, sizeV, devSino, devVoxel);

}
__global__ void yzPlaneFor(const int sizeD[3], const int sizeV[3], float* devSino, const float* devVoxel, GeometryCUDA& geom, const int x, const int n) {
    const int y = blockIdx.x * blockDim.x + threadIdx.x;
    const int z = blockIdx.y * blockDim.y + threadIdx.y;

    const int coord[4] = {x, y, z, n};

    geom.forwardProj(coord, sizeD, sizeV, devSino, devVoxel);
}

__global__ void backwardProj(const int sizeD[3], const int sizeV[3], const float* devSino, float* devVoxel, GeometryCUDA& geom, const int z, const int n) {
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;

    const int coord[4] = {x, y, z, n};

    geom.backwardProj(coord, sizeD, sizeV, devSino, devVoxel);
}


inline void MLEM_CUDA::reconstruct(Volume<float> &sinogram, Volume<float> &voxel, const GeometryCUDA &geom, const int epoch, const int batch, Rotate dir) {
    int sizeV[3] = {voxel.size()[0], voxel.size()[1], voxel.size()[2]};

}