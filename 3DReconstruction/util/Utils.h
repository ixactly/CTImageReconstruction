//
// Created by 森智希 on 2022/06/03.
//

#ifndef INC_3DRECONSTRUCTION_UTILS_H
#define INC_3DRECONSTRUCTION_UTILS_H

#include <array>

enum class Rotate {
    CW,
    CCW
};

using Vec3i = std::array<int, 3>;
using Vec3f = std::array<float, 3>;
using Vec3d = std::array<double, 3>;

struct int3 {
    int sizeX;
    int sizeY;
    int sizeZ;
};

template <typename T>
int sign(T val) {
    return (val > T(0)) - (val < T(0));
}
#endif //INC_3DRECONSTRUCTION_UTILS_H
