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

    std::pair<double, double>
    vox2det(const uint32_t x, const uint32_t y, const uint32_t z, const Vec3i &sizeV, const Vec3i &sizeD,
            double theta) const {
        // impl
        // (x+0.5f, y+0.5f, z+0.5f), source point間の関係からdetのu, vを算出
        Vec3d vecSod = {std::sin(theta) * sod, -std::cos(theta) * sod, 0};
        Vec3d src2cent = {-vecSod[0], -vecSod[1], sizeV[2] * 0.5 * voxSize};
        Vec3d src2vox = {(2 * x - (sizeV[0] - 1)) * 0.5 * voxSize + src2cent[0],
                         (2 * y - (sizeV[1] - 1)) * 0.5 * voxSize + src2cent[1],
                         (2 * z - (sizeV[2] - 1)) * 0.5 * voxSize + src2cent[2]}; // a

        double beta = std::acos(src2vox[0] * src2cent[0] + src2vox[1] * src2cent[1]
                                                           / std::sqrt(
                src2vox[0] * src2vox[0] + src2vox[1] * src2vox[1]) +
                                std::sqrt(src2cent[0] * src2cent[0] + src2cent[1] * src2cent[1]));
        double gamma = std::atan2(src2vox[2], std::sqrt(src2vox[0] * src2vox[0] + src2vox[1] * src2vox[1]));
        double u = std::tan(beta) * sdd / detSize + sizeD[0] * 0.5;
        double v = std::fabs(std::tan(beta)) * std::tan(gamma) * sdd / detSize + sizeD[1] * 0.5; // normalization

        return std::make_pair(u, v);
    }

    bool isHitDetect(const double u, const double v, const Vec3i &sizeD) const {
        // impl
        if (0 < u && u < sizeD[0] && 0 < v && v < sizeD[1])
            return true;
        else
            return false;
    }

private:
    double sdd; // Object-Detector Distance
    double sod; // Source-Object Distance

    double voxSize; // voxel size
    double detSize; // detector size

    Vec3d axisOffset; // point of object rotation center from rotCenter

};

#endif //INC_3DRECONSTRUCTION_GEOMETRY_H
