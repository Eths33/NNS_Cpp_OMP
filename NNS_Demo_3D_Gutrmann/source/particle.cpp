/*
* Author: Gregory Gutmann
* Nearest neighbor search algorithm demo.
*/

#include <particle.hpp>
#include <sort.hpp>
#include <globals.hpp>

void Particle::init(int particleCount, int dimx, int dimy, int dimz) {
	float halfDimx = dimx / 2.0f;
	float halfDimy = dimy / 2.0f;
	float halfDimz = dimz / 2.0f;

	count = particleCount;

	// Place particles at random locations
	for (int i = 0; i < particleCount; i++) {
		float xLoc = ((rand() % (dimx * 10)) / 10.0f) - halfDimx;
		float yLoc = ((rand() % (dimy * 10)) / 10.0f) - halfDimy;
		float zLoc = ((rand() % (dimz * 10)) / 10.0f) - halfDimz;
		locations.push_back(xLoc);
		locations.push_back(yLoc);
		locations.push_back(zLoc);
		neighborCount.push_back(0);
	}

	sortedLoc.resize(locations.size());
	neighborCountN2.resize(neighborCount.size());

	std::vector<std::vector<int>> temp(neighborCount.size());
	neighborList = temp;
	neighborN2List = temp;
}

// All-to-all interaction alogithim O(n^2)
void Particle::countNeighborsN2(int cellLength) {

	int currIdx = 0;

#if PERFORMANCE_TEST && MULTI_THREAD
#pragma omp parallel for
#endif
    for (currIdx = 0; currIdx < count; currIdx++) {
		float3 thisLoc = make_float3(locations[currIdx * 3 + 0], locations[currIdx * 3 + 1], locations[currIdx * 3 + 2]);

        int localCount = 0;

		for (int checkIdx = 0; checkIdx < count; checkIdx++) { 

			if (checkIdx != currIdx) // Dont compute with its self
			{
				float3 checkIdxLoc = make_float3(locations[checkIdx * 3 + 0], locations[checkIdx * 3 + 1], locations[checkIdx * 3 + 2]);
				float3 p2pVec = make_float3(checkIdxLoc.x - thisLoc.x, checkIdxLoc.y - thisLoc.y, checkIdxLoc.z - thisLoc.z);
				float dist = sqrtf((p2pVec.x * p2pVec.x) + (p2pVec.y * p2pVec.y) + (p2pVec.z * p2pVec.z));
				if (dist < (float)cellLength)
				{
					++localCount;
#if !PERFORMANCE_TEST
					neighborN2List[currIdx].push_back(checkIdx);
#endif
				}
			}
		}

        neighborCountN2[currIdx] = localCount;

    }
}

// Using the NNS to run "short" range algorithm
void Particle::countNeighbors(NNS sort) {

	int currIdx = 0;

#if PERFORMANCE_TEST && MULTI_THREAD
#pragma omp parallel for
#endif
	for (currIdx = 0; currIdx < count; currIdx++) {
		float3 thisLoc = make_float3(sortedLoc[currIdx * 3 + 0], sortedLoc[currIdx * 3 + 1], sortedLoc[currIdx * 3 + 2]);
		int thisCell = sort.cellIndexPair[currIdx].cellID;
		int originalIndex = sort.cellIndexPair[currIdx].index;

		int localCount = 0;

		// Check the cells around the particle 
		for (int t = 0; t < 27; t++) {

			// Iterate through the 27 neighbor cells
			int targetCell = (sort.cellDimy * sort.cellDimx * ((t / 9) - 1)) + ((thisCell - sort.cellDimx + (((t % 9) / 3) * sort.cellDimx)) - 1 + ((t) % 3));

			if (targetCell < sort.cellCount - 1) // Excludes the one out-of-bounds cell
			{
				uint32_t startIndex = sort.cellStart[targetCell];
				uint32_t endIndex = sort.cellEnd[targetCell];

				if (startIndex != 0xffffffff)          // cell is not empty
				{
					for (uint32_t checkIdx = startIndex; checkIdx < endIndex; checkIdx++) { // This iterator is going through indexes of the sorted data
						// If data is not sorted will need to use the result of the keyValue sort to get the unsorted (original) index

						if (checkIdx != currIdx) // Dont compute with its self
						{
							float3 checkIdxLoc = make_float3(sortedLoc[checkIdx * 3 + 0], sortedLoc[checkIdx * 3 + 1], sortedLoc[checkIdx * 3 + 2]);

							// Careful with this vectors direction for different interactions
							float3 p2pVec = make_float3(checkIdxLoc.x - thisLoc.x, checkIdxLoc.y - thisLoc.y, checkIdxLoc.z - thisLoc.z);

							float dist = sqrtf((p2pVec.x * p2pVec.x) + (p2pVec.y * p2pVec.y) + (p2pVec.z * p2pVec.z));

							if (dist < (float)sort.cellLength)
							{
								++localCount;
#if !PERFORMANCE_TEST
								// Get checkIdx's unsorted index to access unsorted data
								int checkOrig = sort.cellIndexPair[checkIdx].index;
								neighborList[originalIndex].push_back(checkOrig);
#endif
							}
						}
					}
				}
			}
		}

		neighborCount[originalIndex] = localCount;

	}
}

// Helper
int minimizePrinting(int loop) {
	if (loop > 100) {
		printf("Count is high will only print 10\n");
		loop = 10;
	}
	return loop;
}

void Particle::printLoc(int printCount) {
	int loop = (printCount) ? printCount : count;
	loop = minimizePrinting(loop);
	
	for (int i = 0; i < loop; i++) {
		printf("%.3f %.3f\n", locations[i * 2 + 0], locations[i * 2 + 1]);
	}
}

// Printing NNS results
void Particle::printNeighborCount(int printCount) {
#if PERFORMANCE_TEST
	printNeighborLess(printCount);
#else
	printNeighborMore(printCount);
#endif
}
void Particle::printNeighborLess(int printCount) {
	int loop = (printCount) ? printCount : count;
	loop = minimizePrinting(loop);

	for (int i = 0; i < loop; i++) {
		printf("Particle %d, Neighbor   count %d\n", i, neighborCount[i]);
	}
}
void Particle::printNeighborMore(int printCount) {
	int loop = (printCount) ? printCount : count;
	loop = minimizePrinting(loop);

	for (int i = 0; i < loop; i++) {
		printf("Particle %d, Neighbor count %d    \t { ", i, neighborCount[i]);
		for (int n = 0; n < neighborList[i].size(); n++) {
			printf("%d%s", neighborList[i][n], (n == neighborList[i].size() - 1) ? " }" : ", ");
		}
		printf("\n");
	}
}

// Printing all-to-all results
void Particle::printNeighborN2Count(int printCount) {
#if PERFORMANCE_TEST
	printNeighborN2Less(printCount);
#else
	printNeighborN2More(printCount);
#endif
}
void Particle::printNeighborN2Less(int printCount) {
	int loop = (printCount) ? printCount : count;
	loop = minimizePrinting(loop);

	for (int i = 0; i < loop; i++) {
		printf("Particle %d, NeighborN2 count %d\n", i, neighborCountN2[i]);
	}
}
void Particle::printNeighborN2More(int printCount) {
	int loop = (printCount) ? printCount : count;
	loop = minimizePrinting(loop);

	for (int i = 0; i < loop; i++) {
		printf("Particle %d, NeighborN2 count %d  \t { ", i, neighborCountN2[i]);
		for (int n = 0; n < neighborN2List[i].size(); n++) {
			printf("%d%s", neighborN2List[i][n], (n == neighborN2List[i].size() -1) ? " }" : ", ");
		}
		printf("\n");
	}
}

// Testing: comparing NNS and all-to-all results
void Particle::check() {
	int foundError = 0;
	printf("Checking for differences between NNS and all-to-all\n");

	for (int i = 0; i < count; i++) {
		int firstMissing = 1;
		
		for (int n = 0; n < neighborN2List[i].size(); n++) {
			int found = 0;
			int search = neighborN2List[i][n];
			
			for (int p = 0; p < neighborList[i].size(); p++) {
				if (search == neighborList[i][p]) {
					found = 1;
				}
			}
			
			if (!found) {
				++foundError;
				if (firstMissing) {
					printf("Particle %d\n", i);
				}
				float3 iLoc = make_float3(locations[i * 3 + 0], locations[i * 3 + 1], locations[i * 3 + 2]);
				float3 nLoc = make_float3(locations[search * 3 + 0], locations[search * 3 + 1], locations[search * 3 + 2]);
				float3 p2pVec = make_float3(nLoc.x - iLoc.x, nLoc.y - iLoc.y, nLoc.z - iLoc.z);
				float dist = sqrtf((p2pVec.x * p2pVec.x) + (p2pVec.y * p2pVec.y) + (p2pVec.z * p2pVec.z));
				printf("\tNNS didnt find %d, distance of %.2f\n", search, dist);
			}
		}
	}
	
	if (foundError) {
		printf("\tFound %d errors (some may be mirrors of eachother)\n", foundError);
	}
	else {
		printf("\tSuccess!\n");
	}
}

int Particle::getParticleCount() {
	return count;
}