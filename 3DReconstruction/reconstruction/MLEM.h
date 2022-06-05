//
// Created by 森智希 on 2022/06/03.
//

#ifndef INC_3DRECONSTRUCTION_ELEM_H
#define INC_3DRECONSTRUCTION_ELEM_H

#include "Geometry.h"
#include "Volume.h"

template<typename T>
void MLEM(const Volume<T> &sinogram, Volume<T> &voxel, const Geometry &geom);

#endif //INC_3DRECONSTRUCTION_ELEM_H
