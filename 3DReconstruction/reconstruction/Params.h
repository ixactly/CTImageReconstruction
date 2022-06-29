//
// Created by 森智希 on 2022/06/05.
//

#ifndef INC_3DRECONSTRUCTION_PARAMS_H
#define INC_3DRECONSTRUCTION_PARAMS_H

// params on yuki
/*
inline double SRC_OBJ_DISTANCE = 455.849;
inline double SRC_DETECT_DISTANCE = 1519.739;

inline int NUM_PROJ = 1000;

inline int NUM_DETECT_U = 1024;
inline int NUM_DETECT_V = 1;

inline double DETECTOR_SIZE = 0.4 * 1024 / NUM_DETECT_U;

inline int NUM_VOXEL = 1024;
*/

// params on cfrp

inline double SRC_OBJ_DISTANCE = 1069.0;
inline double SRC_DETECT_DISTANCE = 1450.0;

inline int NUM_PROJ = 360;

inline int NUM_DETECT_U = 1344;
inline int NUM_DETECT_V = 1;

inline double DETECTOR_SIZE = 100.5312 / 1344.0;

inline int NUM_VOXEL = 1344;

#endif //INC_3DRECONSTRUCTION_PARAMS_H
