/*
* Author: Gregory Gutmann
* Nearest neighbor search algorithm demo
*/

#include <sort.hpp>
#include <iostream>
#include <algorithm> // for sort function

KeyValuePair NNS::makeKeyValue(int cell, int idx) {
    KeyValuePair temp;
    temp.cellID = cell;
    temp.index = idx;
    return temp;
}

/// WARNING: A lot of the values around cell assume friendly evenly divisible numbers here
void NNS::init(int count, int dimx, int dimy, int dimz, int cell, int buffer) {
    // Multiply buffer by two to get that amount of buffer on all sides
    simDimx_buffered = dimx + (float)buffer * 2.0f;
    simDimy_buffered = dimy + (float)buffer * 2.0f;
    simDimz_buffered = dimz + (float)buffer * 2.0f;

    //Truncation will likely occure here, be careful
    cellDimx = (int)simDimx_buffered / cell;
    cellDimy = (int)simDimy_buffered / cell;
    cellDimz = (int)simDimz_buffered / cell;

    cellLength = cell;
    bufferSize = buffer;
    cellCount = cellDimx * cellDimy * cellDimz;

    nonBufferCellEstimate = (dimx / cell) * (dimy / cell) * (dimz / cell); // Not used in algorithm

    cellStart.resize(cellCount);
    cellEnd.resize(cellCount);

    particleCount = count;
    cellIndexPair.resize(particleCount);
}

// Utility comparator function to pass to the sort() module
bool sortByCellID(const KeyValuePair& a, const KeyValuePair& b)
{
    return (a.cellID < b.cellID);
}

// In use
void NNS::hash(std::vector<float> &locations) {
    float xShift = simDimx_buffered / 2.0f;
    float yShift = simDimy_buffered / 2.0f;
    float zShift = simDimz_buffered / 2.0f;

    for (int i = 0; i < particleCount; i++){

        int yCube, xCube, zCube;
        xCube = (int)(locations[i * 3 + 0] + xShift) / cellLength;
        yCube = (int)(locations[i * 3 + 1] + yShift) / cellLength;
        zCube = (int)(locations[i * 3 + 2] + zShift) / cellLength;
        int hash = xCube + yCube * cellDimx + zCube * cellDimx * cellDimy;

        if (hash >= cellCount || hash < 0) {
#if defined(DEBUG) | defined(_DEBUG)
            printf("Hash ERROR: HashVal %u - Max HashVal %u   c(%d,%d,%d) l(%f,%f,%f)\n",
                hash, cellCount, xCube, yCube, zCube, locations[i * 3 + 0], locations[i * 3 + 1], locations[i * 3 + 2]);
#endif
            hash = cellCount - 1;		// In calculation this cell is excluded (it is in the outter buffer region)
        }

        cellIndexPair[i].cellID = hash;
        cellIndexPair[i].index = i;
    }
}

// Testing
int NNS::hash(float3 location) {
	float xShift = simDimx_buffered / 2.0f;
	float yShift = simDimy_buffered / 2.0f;
    float zShift = simDimz_buffered / 2.0f;

	int yCube, xCube, zCube;
	xCube = (int)(location.x + xShift) / cellLength;
	yCube = (int)(location.y + yShift) / cellLength;
	zCube = (int)(location.z + zShift) / cellLength;
	int hash = xCube + yCube * cellDimx + zCube * cellDimx * cellDimy;

	if (hash >= cellCount || hash < 0) {
#if defined(DEBUG) | defined(_DEBUG)
        printf("Hash ERROR: HashVal %u - Max HashVal %u   c(%d,%d,%d) l(%f,%f,%f)\n",
            hash, cellCount, xCube, yCube, zCube, location.x, location.y, location.z);
#endif
		hash = cellCount - 1;		// In calculation this cell is excluded (it is in the outter buffer region)
	}

	return hash;
}

void NNS::kvSort() {
    // sort the vector by increasing order of its cell ID
    sort(cellIndexPair.begin(), cellIndexPair.end(), sortByCellID);
}

void NNS::findCellStartEnd() {

    // Mem set to signal cells are empty if not set in this function
    memset(cellStart.data(), 0xffffffff, cellStart.size() * sizeof(uint32_t));

    uint32_t current = 0;
    for (int i = 0; i < particleCount; i++) {
        uint32_t cell = cellIndexPair[i].cellID;
        if (cell != current) { // Found entry in new cell
            cellEnd[current] = i;
            cellStart[cell] = i;
            current = cell;
        }
    }
    cellEnd[current] = particleCount; // Handle last item
}

void NNS::reorder(std::vector<float>& locations, std::vector<float>& sortedLoc) {
    for (int i = 0; i < particleCount; ++i) {
        int originalIndex = cellIndexPair[i].index;

        sortedLoc[i * 3 + 0] = locations[originalIndex * 3 + 0];
        sortedLoc[i * 3 + 1] = locations[originalIndex * 3 + 1];
        sortedLoc[i * 3 + 2] = locations[originalIndex * 3 + 2];
    }
}

// Helper
int minimizePrint(int loop) {
    if (loop > 100) {
        printf("Count is high will only print 10\n");
        loop = 10;
    }
    return loop;
}

// Printing main data strutures used in the NNS
void NNS::printCellIndexPair(int printCount) {
    int loop = (printCount) ? printCount : particleCount;
    loop = minimizePrint(loop);

    for (int i = 0; i < loop; i++) {
        printf("Cell %d, Index %d\n", cellIndexPair[i].cellID, cellIndexPair[i].index);
    }
}
void NNS::printCellStartEnd(int printCount) {
    int loop = (printCount) ? printCount : (int)cellStart.size();
    loop = minimizePrint(loop);

    for (int i = 0; i < loop; i++) {
        if (cellStart[i] == 0xffffffff) {
            printf("Cell %i: Empty\n", i);
        }
        else {
            printf("Cell %i: Start %u, End %u\n", i, cellStart[i], cellEnd[i]);
        }
    }
}

int NNS::getCellCount() {
    return cellCount;
}

int NNS::getNonBuffCellCount() {
    return nonBufferCellEstimate;
}