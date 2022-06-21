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
        voxSize = sod * detSize / sdd;
    }

    std::tuple<double, double, double>
    vox2det(const int x, const int y, const int z, const Vec3i &sizeV, const Vec3i &sizeD,
            double theta) const {
        // impl
        // (x+0.5f, y+0.5f, z+0.5f), source point間の関係からdetのu, vを算出
        // detectorの底辺と光源のz座標は一致していると仮定

        Vec3d vecSod = {std::sin(theta) * sod, -std::cos(theta) * sod, sizeV[2] * 0.5 * voxSize};
        Vec3d src2cent = {-vecSod[0], -vecSod[1], sizeV[2] * 0.5 * voxSize - vecSod[2]};
        Vec3d src2voxel = {(2 * x - sizeV[0] + 1) * 0.5 * voxSize + src2cent[0],
                           (2 * y - sizeV[1] + 1) * 0.5 * voxSize + src2cent[1],
                           (2 * z - sizeV[2] + 1) * 0.5 * voxSize + src2cent[2]}; // a

        double beta = std::acos((src2cent[0] * src2voxel[0] + src2cent[1] * src2voxel[1]) /
                                (std::sqrt(src2cent[0] * src2cent[0] + src2cent[1] * src2cent[1]) *
                                 std::sqrt(src2voxel[0] * src2voxel[0] + src2voxel[1] * src2voxel[1])));
        int signature = sign(src2voxel[0] * src2cent[1] - src2voxel[1] * src2cent[0]); // src2voxel x src2cent
        double gamma = std::atan2(src2voxel[2], std::sqrt(src2voxel[0] * src2voxel[0] + src2voxel[1] * src2voxel[1]));

        double u = std::tan(signature * beta) * sdd / detSize + sizeD[0] * 0.5;
        double v = std::tan(gamma) * sdd / std::cos(beta) / detSize + sizeD[1] * 0.5; // normalization

        return std::make_tuple(u, v, beta);
    }

    bool isHitDetect(const double u, const double v, const Vec3i &sizeD) const {
        // impl
        /*
        if (0.5 < u && u < sizeD[0] - 0.5 && 0.5 < v && v < sizeD[1] - 0.5)
            return true;
        else
            return false;
        */
        if (0.5 < u && u < sizeD[0] - 0.5 && 0 < v && v < sizeD[1])
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
