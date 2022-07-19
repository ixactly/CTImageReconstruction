//
// Created by 森智希 on 2022/06/03.
//

#ifndef INC_3DRECONSTRUCTION_UTILS_H
#define INC_3DRECONSTRUCTION_UTILS_H

#include <array>

template <typename T>
int sign(T val) {
    return (val > T(0)) - (val < T(0));
}

enum class Rotate {
    CW,
    CCW
};

using Vec3i = std::array<int, 3>;
using Vec3f = std::array<float, 3>;
using Vec3d = std::array<double, 3>;

#endif //INC_3DRECONSTRUCTION_UTILS_H
