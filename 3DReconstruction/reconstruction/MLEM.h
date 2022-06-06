//
// Created by 森智希 on 2022/06/03.
//

#ifndef INC_3DRECONSTRUCTION_ELEM_H
#define INC_3DRECONSTRUCTION_ELEM_H

#include "Geometry.h"
#include "Volume.h"
#include "Utils.h"
#include <functional>
#include <utility>
#include <omp.h>

template<typename T>
class MLEM {
public :
    MLEM() = default;

    void reconstruct(const Volume<T> &sinogram, Volume<T> &voxel, const Geometry &geom) {
        Vec3i s = sinogram.size();
        Vec3i vSize = voxel.size();
        Volume<T> volTmp(vSize[0], vSize[1], vSize[2]);
        Volume<T> projTmp(s[0], s[1], s[2]);

        volTmp.forEach([](T value) -> T { return 1.0; });
        projTmp.forEach([](T value) -> T { return 0.0; });

        int nProj = s[2];
        double theta;

        // forward proj
// #pragma omp parallel for
        for (int n = 0; n < nProj; n++) {
            for (int z = 0; z < vSize[2]; z++) {
                for (int y = 0; y < vSize[1]; y++) {
                    for (int x = 0; x < vSize[0]; x++) {
                        theta = 2.0 * M_PI * n / nProj;
                        // forward projection
                        auto [u, v] = geom.vox2det(x, y, z, voxel.size(), sinogram.size(), theta);
                        if (geom.isHitDetect(u, v, sinogram.size())) {
                            lerpScatter(u, v, n, projTmp, voxel(x, y, z)); //need impl
                        }
                    }
                }
            }
        }

        // back proj (y'/y)
        for (int z = 0; z < vSize[2]; z++) {
            for (int y = 0; y < vSize[1]; y++) {
                for (int x = 0; x < vSize[0]; x++) {
                    double c = 0;
                    T voxTmp = 0;
                    for (int n = 0; n < nProj; n++) {
                        theta = 2.0 * M_PI * n / nProj;
                        auto [u, v] = geom.vox2det(x, y, z, voxel.size(), sinogram.size(), theta);
                        if (geom.isHitDetect(u, v, sinogram.size())) {
                            auto [c1, c2, c3, c4] = lerpGather(u, v, n, projTmp, voxTmp);
                            c += c1 + c2 + c3 + c4;
                        }
                    }
                    voxel(x, y, z) = voxTmp * voxel(x, y, z) / c;
                }
            }
        }
    }

    void forwardproj(Volume<T> &sinogram, const Volume<T> &voxel, const Geometry &geom) {
        Vec3i s = sinogram.size();
        Vec3i vSize = voxel.size();
        // Volume<T> volTmp(vSize[0], vSize[1], vSize[2]);
        // Volume<T> projTmp(s[0], s[1], s[2]);

        int nProj = s[2];
        double theta;

        // forward proj
        for (int n = 0; n < nProj; n++) {

            for (int z = 0; z < vSize[2]; z++) {
                for (int y = 0; y < vSize[1]; y++) {
                    for (int x = 0; x < vSize[0]; x++) {
                        theta = 2.0 * M_PI * n / nProj;
                        // forward projection
                        auto [u, v] = geom.vox2det(x, y, z, voxel.size(), sinogram.size(), theta);
                        if (geom.isHitDetect(u, v, sinogram.size())) {
                            lerpScatter(u, v, n, sinogram, voxel(x, y, z)); //need impl
                        }
                    }
                }
            }
        }
    }

    std::tuple<double, double, double, double>
    lerpGather(const double u, const double v, const int n, const Volume<T> &proj,
               T &val) { // normalize u, v on vox2det
        int intU = std::floor(u), intV = std::floor(v);
        double c1 = (1.0 - (u - intU)) * (v - intV), c2 = (u - intU) * (v - intV),
                c3 = (u - intU) * (1.0 - (v - intV)), c4 = (1.0 - (u - intU)) * (1.0 - (v - intV));

        val += c1 * proj(intU, intV + 1, n) + c2 * proj(intU + 1, intV + 1, n) + c3 * proj(intU + 1, intV, n) +
               c4 * proj(intU, intV, n);

        return std::make_tuple(c1, c2, c3, c4);
    }

    std::tuple<double, double, double, double>
    lerpScatter(const double u, const double v, const int n, Volume<T> &proj, const T &val) {
        int intU = std::floor(u), intV = std::floor(v);
        double c1 = (1.0 - (u - intU)) * (v - intV), c2 = (u - intU) * (v - intV),
                c3 = (u - intU) * (1.0 - (v - intV)), c4 = (1.0 - (u - intU)) * (1.0 - (v - intV));

        proj(intU, intV + 1, n) += c1 * val;
        proj(intU + 1, intV + 1, n) += c2 * val;
        proj(intU + 1, intV, n) += c3 * val;
        proj(intU, intV, n) += c4 * val;

        return std::make_tuple(c1, c2, c3, c4);
    }
};


#endif //INC_3DRECONSTRUCTION_ELEM_H
