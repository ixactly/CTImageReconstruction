//
// Created by 森智希 on 2022/06/01.
//

#ifndef INC_3DRECONSTRUCTION_VOLUME_H
#define INC_3DRECONSTRUCTION_VOLUME_H

#include <memory>
#include <array>
#include <string>
#include <opencv2/opencv.hpp>

template<typename T>
class Volume {
public :
    Volume() = default;

    explicit Volume(uint32_t sizeX, uint32_t sizeY, uint32_t sizeZ)
            : sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ) {}

    explicit Volume(std::string filename, uint32_t sizeX, uint32_t sizeY, uint32_t sizeZ) {
        // implement
    }

    Volume(const Volume &v)
            : sizeX(v.sizeX), sizeY(v.sizeY), sizeZ(v.sizeZ) {
        const uint32_t size = v.sizeX * v.sizeY * v.sizeZ;
        std::memcpy(data.get(), v.data.get(), size * sizeof(uint32_t));
    }

    Volume &operator=(const Volume &v) {
        sizeX = v.sizeX;
        sizeY = v.sizeY;
        sizeZ = v.sizeZ;
        const uint32_t size = v.sizeX * v.sizeY * v.sizeZ;
        std::memcpy(data.get(), v.data.get(), size * sizeof(uint32_t));
    }

    Volume(Volume &&v) : sizeX(v.sizeX), sizeY(v.sizeY), sizeZ(v.sizeZ) {
        data = std::move(v.data);
    }

    Volume &operator=(Volume &&v) {
        sizeX = v.sizeX;
        sizeY = v.sizeY;
        sizeZ = v.sizeZ;
        data = std::move(v.data);

        return *this;
    }

    ~Volume() = default;

    // show the slice of center
    void slice() {

    }


    void show() {

    }

private :
    uint32_t sizeX;
    uint32_t sizeY;
    uint32_t sizeZ;
    std::unique_ptr<T[]> data = nullptr;
};

#endif //INC_3DRECONSTRUCTION_VOLUME_H
