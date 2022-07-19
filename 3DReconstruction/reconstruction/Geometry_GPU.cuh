//
// Created by tomokimori on 22/07/18.
//

#ifndef INC_3DRECONSTRUCTION_GEOMETRY_GPU_CUH
#define INC_3DRECONSTRUCTION_GEOMETRY_GPU_CUH
#include "cuda_runtime.h"
#include "cuda.h"
#include <cmath>

template <typename T>
__host__ __device__ int dsign(T val) {
    return (val > T(0)) - (val < T(0));
}

class GeometryCUDA {

public:
    GeometryCUDA(float sdd, float sod, float detSize) : sdd(sdd), sod(sod), detSize(detSize) {
        voxSize = sod * detSize / sdd;
    }

    __device__ void
    forwardProj(const int coord[4], const int sizeD[3], const int sizeV[3], float *devSino, const float* devVoxel) const {

        // sourceとvoxel座標間の関係からdetのu, vを算出
        // detectorの中心 と 再構成領域の中心 と 光源 のz座標は一致していると仮定
        const int n = coord[3];
        const float theta = 2.0f * M_PI * n / sizeD[2];

        float offset[3] = {0.0, 0.0, 0.0};
        float vecSod[3] = {std::sin(theta) * sod + offset[0], -std::cos(theta) * sod + offset[1], 0};

        // Source to voxel center
        float src2cent[3] = {-vecSod[0], -vecSod[1], -vecSod[2]};
        // Source to voxel
        float src2voxel[3] = {(2 * coord[0] - sizeV[0] + 1) * 0.5f * voxSize + src2cent[0],
                           (2 * coord[1] - sizeV[1] + 1) * 0.5f * voxSize + src2cent[1],
                           (2 * coord[2] - sizeV[2] + 1) * 0.5f * voxSize + src2cent[2]};

        const float beta = std::acos((src2cent[0] * src2voxel[0] + src2cent[1] * src2voxel[1]) /
                                      (std::sqrt(src2cent[0] * src2cent[0] + src2cent[1] * src2cent[1]) *
                                       std::sqrt(src2voxel[0] * src2voxel[0] + src2voxel[1] * src2voxel[1])));
        const float gamma = std::atan2(src2voxel[2], std::sqrt(src2voxel[0]*src2voxel[0]+src2voxel[1]*src2voxel[1]));

        const int signU = dsign(src2voxel[0] * src2cent[1] - src2voxel[1] * src2cent[0]);

        // src2voxel x src2cent
        // 光線がhitするdetector平面座標の算出(detectorSizeで除算して、正規化済み)
        float u = std::tan(signU * beta) * sdd / detSize + sizeD[0] * 0.5f;
        float v = std::tan(gamma) * sdd / std::cos(beta) / detSize + sizeD[1] * 0.5f; // normalization

        if (!(0.5 < u && u < sizeD[0] - 0.5 && 0.5 < v && v < sizeD[1] - 0.5))
            return;

        float u_tmp = u - 0.5f, v_tmp = v - 0.5f;
        int intU = std::floor(u_tmp), intV = std::floor(v_tmp);
        float c1 = (1.0f - (u_tmp - intU)) * (v_tmp - intV), c2 = (u_tmp - intU) * (v_tmp - intV),
                c3 = (u_tmp - intU) * (1.0f - (v_tmp - intV)), c4 =
                (1.0f - (u_tmp - intU)) * (1.0f - (v_tmp - intV));

        const unsigned int idxVoxel = coord[0] + sizeV[0] * coord[1] + sizeV[0] * sizeV[1] * coord[2];

        devSino[intU + sizeD[0] * (intV+1) + sizeD[0] * sizeD[1] * n] += c1 * devVoxel[idxVoxel];
        devSino[(intU+1) + sizeD[0] * (intV+1) + sizeD[0] * sizeD[1] * n] += c2 * devVoxel[idxVoxel];
        devSino[(intU+1) + sizeD[0] * intV + sizeD[0] * sizeD[1] * n] += c3 * devVoxel[idxVoxel];
        devSino[intU + sizeD[0] * intV + sizeD[0] * sizeD[1] * n] += c4 * devVoxel[idxVoxel];
    }

    __device__ void backwardProj(const int coord[4], const int sizeD[3], const int sizeV[3], const float *devSino, float* devVoxel) {
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

private:
    float sdd; // Object-Detector Distance
    float sod; // Source-Object Distance

    float voxSize; // voxel size
    float detSize; // detector size

    // Vec3d axisOffset; // point of object rotation center from rotCenter

};
#endif //INC_3DRECONSTRUCTION_GEOMETRY_GPU_CUH
