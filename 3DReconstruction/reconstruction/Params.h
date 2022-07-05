//
// Created by 森智希 on 2022/06/05.
//

#ifndef INC_3DRECONSTRUCTION_PARAMS_H
#define INC_3DRECONSTRUCTION_PARAMS_H

// params on yuki
/*
inline constexpr double SRC_OBJ_DISTANCE = 455.849;
inline constexpr double SRC_DETECT_DISTANCE = 1519.739;

inline constexpr int NUM_PROJ = 1000;

inline constexpr int NUM_DETECT_U = 1024;
inline constexpr int NUM_DETECT_V = 1024;

inline constexpr double DETECTOR_SIZE = 0.4 * 1024 / NUM_DETECT_U;

inline constexpr int NUM_VOXEL = 1024;
*/

// params on cfrp
/*
inline constexpr double SRC_OBJ_DISTANCE = 1069.0;
inline constexpr double SRC_DETECT_DISTANCE = 1450.0;

inline constexpr int NUM_PROJ = 360;

inline constexpr int NUM_DETECT_U = 1344;
inline constexpr int NUM_DETECT_V = 1;

inline constexpr double DETECTOR_SIZE = 100.5312 / 1344.0;

inline constexpr int NUM_VOXEL = 1344;
*/

inline constexpr double SRC_OBJ_DISTANCE = 500.0;
inline constexpr double SRC_DETECT_DISTANCE = 1500.0;

inline constexpr int NUM_PROJ = 100;

inline constexpr int NUM_DETECT_U = 100;
inline constexpr int NUM_DETECT_V = 100;

inline constexpr double DETECTOR_SIZE = 0.1;

inline constexpr int NUM_VOXEL = 100;
#endif //INC_3DRECONSTRUCTION_PARAMS_H
