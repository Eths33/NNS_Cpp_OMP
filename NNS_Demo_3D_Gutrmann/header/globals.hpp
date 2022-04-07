/*
* Author: Gregory Gutmann
* Nearest neighbor search algorithm demo.
*/

#ifndef HELPER_H
#define HELPER_H

#define PERFORMANCE_TEST 0
#define MULTI_THREAD 1

#if !PERFORMANCE_TEST
    #define X_DIM 10
    #define Y_DIM 10
    #define Z_DIM 5
    #define CELL_SIZE 5
    #define GRID_BUFFER 5 // two times the cell size is a good starting point
    #define PARTICLE_COUNT 10
#elif (PERFORMANCE_TEST && !MULTI_THREAD)
    #define X_DIM 40
    #define Y_DIM 40
    #define Z_DIM 30
    #define CELL_SIZE 5
    #define GRID_BUFFER 10 // two times the cell size is a good starting point
    #define PARTICLE_COUNT 800
#else
    #define X_DIM 60
    #define Y_DIM 60
    #define Z_DIM 60
    #define CELL_SIZE 5
    #define GRID_BUFFER 10 // two times the cell size is a good starting point
    #define PARTICLE_COUNT 3600
#endif

struct float3 {
    float x;
    float y;
    float z;
};

float3 make_float3(float a, float b, float c);

#endif // HELPER_H