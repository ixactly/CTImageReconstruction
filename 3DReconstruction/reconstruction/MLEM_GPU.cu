//
// Created by tomokimori on 22/07/16.
//

#include "MLEM_GPU.cuh"
#include "cuda.h"
#include "cuda_runtime.h"
#include "Geometry_GPU.cuh"
#include "../util/Pbar.h"
#include <random>

__device__ void GeometryCUDA::forwardProj(const int coord[4], const int sizeD[3], const int sizeV[3], float *devSino, const float* devVoxel) const {

    // sourceとvoxel座標間の関係からdetのu, vを算出
    // detectorの中心 と 再構成領域の中心 と 光源 のz座標は一致していると仮定
    const int n = coord[3];
    const float theta = 2.0f * M_PI * n / sizeD[2];

    float offset[3] = {0.0, 0.0, 0.0};
    float vecSod[3] = {sinf(theta) * sod + offset[0], -cosf(theta) * sod + offset[1], 0};

    // Source to voxel center
    float src2cent[3] = {-vecSod[0], -vecSod[1], -vecSod[2]};
    // Source to voxel
    float src2voxel[3] = {(2.0f * coord[0] - sizeV[0] + 1) * 0.5f * voxSize + src2cent[0],
                          (2.0f * coord[1] - sizeV[1] + 1) * 0.5f * voxSize + src2cent[1],
                          (2.0f * coord[2] - sizeV[2] + 1) * 0.5f * voxSize + src2cent[2]};

    const float beta = acos((src2cent[0] * src2voxel[0] + src2cent[1] * src2voxel[1]) /
                            (sqrt(src2cent[0] * src2cent[0] + src2cent[1] * src2cent[1]) *
                             sqrt(src2voxel[0] * src2voxel[0] + src2voxel[1] * src2voxel[1])));
    const float gamma = atan2(src2voxel[2], sqrt(src2voxel[0]*src2voxel[0]+src2voxel[1]*src2voxel[1]));

    const int signU = dsign(src2voxel[0] * src2cent[1] - src2voxel[1] * src2cent[0]);

    // src2voxel x src2cent
    // 光線がhitするdetector平面座標の算出(detectorSizeで除算して、正規化済み)
    float u = tanf(signU * beta) * sdd / detSize + sizeD[0] * 0.5f;
    float v = tanf(gamma) * sdd / cosf(beta) / detSize + sizeD[1] * 0.5f; // normalization

    if (!(0.5 < u && u < sizeD[0] - 0.5 && 0.5 < v && v < sizeD[1] - 0.5))
        return;

    float u_tmp = u - 0.5f, v_tmp = v - 0.5f;
    int intU = floor(u_tmp), intV = floor(v_tmp);
    float c1 = (1.0f - (u_tmp - intU)) * (v_tmp - intV), c2 = (u_tmp - intU) * (v_tmp - intV),
            c3 = (u_tmp - intU) * (1.0f - (v_tmp - intV)), c4 =
            (1.0f - (u_tmp - intU)) * (1.0f - (v_tmp - intV));

    const unsigned int idxVoxel = coord[0] + sizeV[0] * coord[1] + sizeV[0] * sizeV[1] * coord[2];

    atomicAdd(&devSino[intU + sizeD[0] * (intV+1) + sizeD[0] * sizeD[1] * n], c1 * devVoxel[idxVoxel]);
    atomicAdd(&devSino[(intU+1) + sizeD[0] * (intV+1) + sizeD[0] * sizeD[1] * n], c2 * devVoxel[idxVoxel]);
    atomicAdd(&devSino[(intU+1) + sizeD[0] * intV + sizeD[0] * sizeD[1] * n], c3 * devVoxel[idxVoxel]);
    atomicAdd(&devSino[intU + sizeD[0] * intV + sizeD[0] * sizeD[1] * n], c4 * devVoxel[idxVoxel]);

    /*
    devSino[intU + sizeD[0] * (intV+1) + sizeD[0] * sizeD[1] * n] += c1 * devVoxel[idxVoxel];
    devSino[(intU+1) + sizeD[0] * (intV+1) + sizeD[0] * sizeD[1] * n] += c2 * devVoxel[idxVoxel];
    devSino[(intU+1) + sizeD[0] * intV + sizeD[0] * sizeD[1] * n] += c3 * devVoxel[idxVoxel];
    devSino[intU + sizeD[0] * intV + sizeD[0] * sizeD[1] * n] += c4 * devVoxel[idxVoxel];
    */
}
__device__ void GeometryCUDA::backwardProj(const int coord[4], const int sizeD[3], const int sizeV[3], const float *devSino, float* devVoxel) const {
    // sourceとvoxel座標間の関係からdetのu, vを算出
    // detectorの中心 と 再構成領域の中心 と 光源 のz座標は一致していると仮定
    const int n = coord[3];
    const float theta = 2 * M_PI * n / sizeD[2];

    float offset[3] = {0.0, 0.0, 0.0};
    float vecSod[3] = {std::sin(theta) * sod + offset[0], -std::cos(theta) * sod + offset[1], 0};

    // Source to voxel center
    float src2cent[3] = {-vecSod[0], -vecSod[1], -vecSod[2]};
    // Source to voxel
    float src2voxel[3] = {(2.0f * coord[0] - sizeV[0] + 1) * 0.5f * voxSize + src2cent[0],
                          (2.0f * coord[1] - sizeV[1] + 1) * 0.5f * voxSize + src2cent[1],
                          (2.0f * coord[2] - sizeV[2] + 1) * 0.5f * voxSize + src2cent[2]};

    const float beta = std::acos((src2cent[0] * src2voxel[0] + src2cent[1] * src2voxel[1]) /
                                 (std::sqrt(src2cent[0] * src2cent[0] + src2cent[1] * src2cent[1]) *
                                  std::sqrt(src2voxel[0] * src2voxel[0] + src2voxel[1] * src2voxel[1])));
    const float gamma = std::atan2(src2voxel[2], std::sqrt(src2voxel[0]*src2voxel[0]+src2voxel[1]*src2voxel[1]));

    const int signU = dsign(src2voxel[0] * src2cent[1] - src2voxel[1] * src2cent[0]);

    // src2voxel x src2cent
    // 光線がhitするdetector平面座標の算出(detectorSizeで除算して、正規化済み)
    float u = std::tan(signU * beta) * sdd / detSize + sizeD[0] * 0.5;
    float v = std::tan(gamma) * sdd / std::cos(beta) / detSize + sizeD[1] * 0.5; // normalization

    if (!(0.5 < u && u < sizeD[0] - 0.5 && 0.5 < v && v < sizeD[1] - 0.5))
        return;

    float u_tmp = u - 0.5, v_tmp = v - 0.5;
    int intU = std::floor(u_tmp), intV = std::floor(v_tmp);
    float c1 = (1.0 - (u_tmp - intU)) * (v_tmp - intV), c2 = (u_tmp - intU) * (v_tmp - intV),
            c3 = (u_tmp - intU) * (1.0 - (v_tmp - intV)), c4 =
            (1.0 - (u_tmp - intU)) * (1.0 - (v_tmp - intV));

    const unsigned int idxVoxel = coord[0] + sizeV[0] * coord[1] + sizeV[0] * sizeV[1] * coord[2];
    devVoxel[idxVoxel] += c1 * devSino[intU + sizeD[0] * (intV+1) + sizeD[0] * sizeD[1] * n]
                          + c2 * devSino[(intU+1) + sizeD[0] * (intV+1) + sizeD[0] * sizeD[1] * n]
                          + c3 * devSino[(intU+1) + sizeD[0] * intV + sizeD[0] * sizeD[1] * n]
                          + c4 * devSino[intU + sizeD[0] * intV + sizeD[0] * sizeD[1] * n];
}

__global__ void //add atomic
xzPlaneFor(const int sizeD[3], const int sizeV[3], float *devSino, const float *devVoxel, GeometryCUDA *geom,
           const int y, const int n) {
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int z = blockIdx.y * blockDim.y + threadIdx.y;

    const int coord[4] = {x, y, z, n};
    devSino[x+sizeD[0]*y+sizeD[1]*sizeD[0]*z] = 1.0f;
    geom->forwardProj(coord, sizeD, sizeV, devSino, devVoxel);

}

__global__ void
yzPlaneFor(const int sizeD[3], const int sizeV[3], float *devSino, const float *devVoxel, const GeometryCUDA &geom,
           const int x, const int n) {
    const int y = blockIdx.x * blockDim.x + threadIdx.x;
    const int z = blockIdx.y * blockDim.y + threadIdx.y;

    const int coord[4] = {x, y, z, n};

    geom.forwardProj(coord, sizeD, sizeV, devSino, devVoxel);
}

__global__ void
backwardProj(const int sizeD[3], const int sizeV[3], const float *devSino, float *devVoxel, const GeometryCUDA &geom,
             const int z, const int n) {
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;

    const int coord[4] = {x, y, z, n};

    geom.backwardProj(coord, sizeD, sizeV, devSino, devVoxel);
}


void
MLEM_CUDA::reconstruct(Volume<float> &sinogram, Volume<float> &voxel, const GeometryCUDA* geom, const int epoch,
                       const int batch, Rotate dir) {
    int sizeV[3] = {voxel.size()[0], voxel.size()[1], voxel.size()[2]};
    int sizeD[3] = {sinogram.size()[0], sinogram.size()[1], sinogram.size()[2]};
    int nProj = sizeD[2];

    float *devSino, *devVoxel;
    GeometryCUDA *devGeom;

    cudaMalloc(&devSino, sizeof(float) * sizeD[0] * sizeD[1] * sizeD[2]);
    cudaMalloc(&devVoxel, sizeof(float) * sizeV[0] * sizeV[1] * sizeV[2]);
    cudaMalloc(&devGeom, sizeof(GeometryCUDA));

    cudaMemcpy(devSino, sinogram.getPtr(), sizeof(float) * sizeD[0] * sizeD[1] * sizeD[2], cudaMemcpyHostToDevice);
    cudaMemcpy(devVoxel, voxel.getPtr(), sizeof(float) * sizeV[0] * sizeV[1] * sizeV[2], cudaMemcpyHostToDevice);
    cudaMemcpy(devGeom, geom, sizeof(GeometryCUDA), cudaMemcpyHostToDevice);

    const int blockSize = 8;
    dim3 block(blockSize, blockSize, 1);
    dim3 grid((sizeV[0] + blockSize - 1) / blockSize, (sizeV[0] + blockSize - 1) / blockSize, 1);

    // forward, divide, backward proj
    int subsetSize = (nProj + batch - 1) / batch;
    std::vector<int> subsetOrder(batch);
    for (int i = 0; i < batch; i++) {
        subsetOrder[i] = i;
    }

    std::mt19937_64 get_rand_mt; // fixed seed
    std::shuffle(subsetOrder.begin(), subsetOrder.end(), get_rand_mt);

    // progress bar
    progressbar pbar(epoch * nProj);

    // main routine

    for (int ep = 0; ep < epoch; ep++) {
        for (int &sub: subsetOrder) {

            for (int subOrder = 0; subOrder < subsetSize; subOrder++) {
                pbar.update();
                int n = (sub + batch * subOrder) % nProj;
                for (int y = 0; y < sizeV[1]; y++) {
                    xzPlaneFor<<<grid, block>>>(sizeD, sizeV, devSino, devVoxel, devGeom, y, n);
                    cudaDeviceSynchronize();
                    cudaGetLastError();
                }
            }
        }
    }

    cudaMemcpy(voxel.getPtr(), devVoxel, sizeof(float) * sizeV[0] * sizeV[1] * sizeV[2], cudaMemcpyDeviceToHost);
    cudaMemcpy(sinogram.getPtr(), devSino, sizeof(float) * sizeD[0] * sizeD[1] * sizeD[2], cudaMemcpyDeviceToHost);

    cudaFree(devSino);
    cudaFree(devVoxel);
    cudaFree(devGeom);

}