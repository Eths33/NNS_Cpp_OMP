/*
* Author: Gregory Gutmann
* Nearest neighbor search algorithm demo
*/

#ifndef SORT_H
#define SORT_H

#include <globals.hpp>
#include <vector>

struct KeyValuePair {
    int cellID;    // Grid cell
    int index;     // Particle index
};

class NNS {
private:
    // Depending on use case make more things private and use getters and setters

public:    
    int cellLength;
    int bufferSize; // as a buffer and to handle truncation from dividing by cell size problem

    int cellCount;
    int nonBufferCellEstimate;

    int cellDimx;
    int cellDimy;
    int cellDimz;

    float simDimx_buffered;
    float simDimy_buffered;
    float simDimz_buffered;

    std::vector<uint32_t> cellStart;
    std::vector<uint32_t> cellEnd;

    int particleCount;

    KeyValuePair makeKeyValue(int cell, int idx);

    std::vector<KeyValuePair> cellIndexPair;

    void init(int count, int dimx, int dimy, int dimz, int cell, int buffer);

    void hash(std::vector<float>& locations);
    int hash(float3 location);
    void kvSort();
    void findCellStartEnd();
    void reorder(std::vector<float>& locations, std::vector<float>& sortedLoc);

    // Printing main data strutures used in the NNS
    void printCellIndexPair(int printCount = 0);
    void printCellStartEnd(int printCount = 0);

    int getCellCount();
    int getNonBuffCellCount();
};

#endif // SORT_H