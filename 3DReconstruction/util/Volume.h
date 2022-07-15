//
// Created by 森智希 on 2022/06/01.
//

#ifndef INC_3DRECONSTRUCTION_VOLUME_H
#define INC_3DRECONSTRUCTION_VOLUME_H

#include <memory>
#include <array>
#include <string>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <functional>
#include "Utils.h"

template<typename T>
class Volume {
public :
    Volume() = default;

    explicit Volume(int sizeX, int sizeY, int sizeZ)
            : sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ) {
        data = std::make_unique<T[]>(sizeX * sizeY * sizeZ);
    }

    explicit Volume(std::string &filename, int sizeX, int sizeY, int sizeZ)
            : sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ) {
        // implement
        load(filename, sizeX, sizeY, sizeZ);
    }

    Volume(const Volume &v)
            : sizeX(v.sizeX), sizeY(v.sizeY), sizeZ(v.sizeZ) {
        const int size = v.sizeX * v.sizeY * v.sizeZ;
        std::memcpy(data.get(), v.data.get(), size * sizeof(T));
    }

    Volume &operator=(const Volume &v) {
        sizeX = v.sizeX, sizeY = v.sizeY, sizeZ = v.sizeZ;
        const int size = v.sizeX * v.sizeY * v.sizeZ;
        std::memcpy(data.get(), v.data.get(), size * sizeof(T));

        return *this;
    }

    Volume(Volume &&v) noexcept: sizeX(v.sizeX), sizeY(v.sizeY), sizeZ(v.sizeZ) {
        v.sizeX = 0, v.sizeY = 0, v.sizeZ = 0;
        data = std::move(v.data);
    }

    Volume &operator=(Volume &&v) noexcept {
        sizeX = v.sizeX, sizeY = v.sizeY, sizeZ = v.sizeZ;
        v.sizeX = 0, v.sizeY = 0, v.sizeZ = 0;
        data = std::move(v.data);

        return *this;
    }

    ~Volume() = default;

    // ref data (mutable)
    T &operator()(int x, int y, int z) {
        return data[z * (sizeX * sizeY) + y * (sizeX) + x];
    }

    T operator()(int x, int y, int z) const {
        return data[z * (sizeX * sizeY) + y * (sizeX) + x];
    }

    // show the slice of center
    void show(const int slice) { // opencv and unique ptr(need use shared ptr?)
        // axis決めるのはメモリの並び的にだるいっす.
        cv::Mat xyPlane(sizeX, sizeY, cv::DataType<T>::type, data.get() + slice * (sizeX * sizeY));

        cv::imshow("slice", xyPlane);
        cv::waitKey(0);
    }

    void load(const std::string &filename, const int x, const int y, const int z) {
        // impl
        sizeX = x, sizeY = y, sizeZ = z;
        const int size = x * y * z;
        data.reset();
        data = std::make_unique<T[]>(size);
        std::ifstream ifile(filename, std::ios::binary);

        ifile.read(reinterpret_cast<char *>(data.get()), sizeof(T) * size);
    }

    void save(const std::string &filename) {
        const int size = sizeX * sizeY * sizeZ;
        std::ofstream ofs(filename, std::ios::binary);
        if (!ofs) {
            std::cout << "file not opened" << std::endl;
            return;
        }
        ofs.write(reinterpret_cast<char *>(data.get()), sizeof(T) * size);
    }

    void transpose() {
        // impl axis swap
        // use std::swap to data
    }

    void forEach(const std::function<T(T)> &f) {
        for (int z = 0; z < sizeZ; z++) {
            for (int y = 0; y < sizeY; y++) {
                for (int x = 0; x < sizeX; x++) {
                    (*this)(x, y, z) = f((*this)(x, y, z));
                }
            }
        }
    }

    Vec3i size() const {
        Vec3i vec = {static_cast<int>(sizeX), static_cast<int>(sizeY), static_cast<int>(sizeZ)};
        return vec;
    }

private :
    int sizeX, sizeY, sizeZ;
    std::unique_ptr<T[]> data = nullptr;
};

#endif //INC_3DRECONSTRUCTION_VOLUME_H
