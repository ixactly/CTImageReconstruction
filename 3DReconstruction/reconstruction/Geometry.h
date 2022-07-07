//
// Created by 森智希 on 2022/06/02.
//

#ifndef INC_3DRECONSTRUCTION_GEOMETRY_H
#define INC_3DRECONSTRUCTION_GEOMETRY_H

#include <string>
#include <array>
#include <cmath>
#include "Utils.h"

class Geometry {

public:
    Geometry(double sdd, double sod, double detSize) : sdd(sdd), sod(sod), detSize(detSize) {
        /* //not impl axis offset
        Vec3d coordSo = {-axisOffset[0], -(axisOffset[1] - sod), -axisOffset[2]};
        // offset(0, sod), (-axisOffset[0], sod - axisOffset[1])
        double sodCoeff = sod * (sod - axisOffset[1]) / sod *
                          std::sqrt(std::pow(axisOffset[0], 2) + std::pow((sod - axisOffset[1]), 2));
        sod = sodCoeff * std::sqrt(std::pow(axisOffset[0], 2) + std::pow((sod - axisOffset[1]), 2));
        voxSize = sod * detSize / sdd;
        */
        voxSize = sod * detSize * 8.0 / sdd;
    }

    [[nodiscard]] std::pair<double, double>
    vox2det(const int x, const int y, const int z, const int n, const Vec3i &sizeD, const Vec3i &sizeV
    ) const {
        // impl

        // sourceとvoxel座標間の関係からdetのu, vを算出
        // detectorの中心 と 再構成領域の中心 と 光源 のz座標は一致していると仮定
        const double theta = 2 * M_PI * n / sizeD[2];
        Vec3d offset = {0.0, 0.0, 0.0};
        Vec3d vecSod = {std::sin(theta) * sod + offset[0], -std::cos(theta) * sod + offset[1], 0};

        // Source to voxel center
        Vec3d src2cent = {-vecSod[0], -vecSod[1], -vecSod[2]};
        // Source to voxel
        Vec3d src2voxel = {(2 * x - sizeV[0] + 1) * 0.5 * voxSize + src2cent[0],
                           (2 * y - sizeV[1] + 1) * 0.5 * voxSize + src2cent[1],
                           (2 * z - sizeV[2] + 1) * 0.5 * voxSize + src2cent[2]}; // a

        const double beta = std::acos((src2cent[0] * src2voxel[0] + src2cent[1] * src2voxel[1]) /
                                      (std::sqrt(src2cent[0] * src2cent[0] + src2cent[1] * src2cent[1]) *
                                       std::sqrt(src2voxel[0] * src2voxel[0] + src2voxel[1] * src2voxel[1])));
        const double gamma = std::atan2(src2voxel[2], std::sqrt(src2voxel[0]*src2voxel[0]+src2voxel[1]*src2voxel[1]));

        const int signU = sign(src2voxel[0] * src2cent[1] - src2voxel[1] * src2cent[0]);

        // src2voxel x src2cent
        // 光線がhitするdetector平面座標の算出(detectorSizeで除算して、正規化済み)
        double u = std::tan(signU * beta) * sdd / detSize + sizeD[0] * 0.5;
        double v = std::tan(gamma) * sdd / std::cos(beta) / detSize + sizeD[1] * 0.5; // normalization

        return std::make_pair(u, v);
    }

    bool isHitDetect(const double u, const double v, const Vec3i &sizeD) const {
        // impl
        if (0.5 < u && u < sizeD[0] - 0.5 && 0.5 < v && v < sizeD[1] - 0.5)
            return true;
        else
            return false;
    }

private:
    double sdd; // Object-Detector Distance
    double sod; // Source-Object Distance

    double voxSize; // voxel size
    double detSize; // detector size

    // Vec3d axisOffset; // point of object rotation center from rotCenter

};

#endif //INC_3DRECONSTRUCTION_GEOMETRY_H
