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

    __device__ void forwardProj(const int coord[4], const int sizeD[3], const int sizeV[3], float *devSino, const float* devVoxel) const;
    __device__ void backwardProj(const int coord[4], const int sizeD[3], const int sizeV[3], const float *devSino, float* devVoxel) const;

private:
    float sdd; // Object-Detector Distance
    float sod; // Source-Object Distance

    float voxSize; // voxel size
    float detSize; // detector size

    // Vec3d axisOffset; // point of object rotation center from rotCenter

};
#endif //INC_3DRECONSTRUCTION_GEOMETRY_GPU_CUH
