//
// Created by 森智希 on 2022/06/05.
//

#ifndef INC_3DRECONSTRUCTION_PARAMS_H
#define INC_3DRECONSTRUCTION_PARAMS_H

inline double SRC_OBJ_DISTANCE = 500;
inline double SRC_DETECT_DISTANCE = 1000;

inline int NUM_DETECT_U = 500;
inline int NUM_DETECT_V = 1;

inline double DETECTOR_SIZE = 0.2 * 500 / NUM_DETECT_U;
inline int NUM_PROJ = 500;
inline int NUM_VOXEL = 500;

#endif //INC_3DRECONSTRUCTION_PARAMS_H
