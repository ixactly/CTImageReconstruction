//
// Created by 森智希 on 2022/06/01.
//

#include <iostream>
#include <Volume.h>

int main() {
    Volume<float> sinogram;
    Volume<float> ct;

    sinogram.load("../volume_bin/sphere-tori-float-500x500x500.raw", 500, 500, 500);
    sinogram.show(250);


}