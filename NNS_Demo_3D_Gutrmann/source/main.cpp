/*
* Author: Gregory Gutmann
* Nearest neighbor search algorithm demo using a 3D space
* Don't print or keep track of neighbors in a list during simulaion
*/

#include <globals.hpp>
#include <sort.hpp>
#include <particle.hpp>
#include <time.h>
#include <omp.h>

int main() {
#if PERFORMANCE_TEST && MULTI_THREAD
	printf("omp_get_max_threads() = %d\n", omp_get_max_threads());
	int numThreads = omp_get_max_threads() / 2;
	omp_set_num_threads(numThreads);
	printf("omp_get_max_threads() = %d\n\n", omp_get_max_threads());
#endif

	// Feel free to change these values to test
	int xDimension = X_DIM;
	int yDimension = Y_DIM;
	int zDimension = Z_DIM;
	int cellSize = CELL_SIZE;
	int gridBuffer = GRID_BUFFER;
	int particleCount = PARTICLE_COUNT;


	NNS sortObject;
	Particle partObject;

	sortObject.init(particleCount, xDimension, yDimension, zDimension, cellSize, gridBuffer);
	partObject.init(particleCount, xDimension, yDimension, zDimension);

	// --- Simulation loop starts here ----------------------------------------------------
	sortObject.hash(partObject.locations);

	sortObject.printCellIndexPair(); printf("\n\n");
	sortObject.kvSort();
	sortObject.printCellIndexPair(); printf("\n\n");

	sortObject.findCellStartEnd();
	sortObject.printCellStartEnd(); printf("\n\n");

	// Improves memory access pattern, also the algorithm functions in sorted order
	sortObject.reorder(partObject.locations, partObject.sortedLoc);


	partObject.countNeighbors(sortObject);
	partObject.printNeighborCount(); printf("\n\n");
	// --- Simulation loop ends here ------------------------------------------------------


	
#if !PERFORMANCE_TEST
	// Testing (Debug)
	partObject.countNeighborsN2(cellSize);
	partObject.printNeighborN2Count(); printf("\n\n");
	partObject.check(); printf("\n\n");
#endif

#if PERFORMANCE_TEST
	clock_t t;
	float nnsTime, ataTime;
	// Running 1000x to get a larger time for the timer, and to get a more repeatable performance number
	printf("Running 1000 iterations of NNS and all-to-all\n\n");
	
	printf("Simulation space is about: %d x %d x %d\n", xDimension, yDimension, zDimension);
	printf("Particle count:            %d\n\n", particleCount);
	printf("Average particles per non buffer cell: %.2f (2.0 is not sparse) \n", partObject.getParticleCount() / (float)sortObject.getNonBuffCellCount());
	printf("At greater densities the NNS performance gain will decrease\n\n");

	// NNS
	{
		t = clock();
		for (int i = 0; i < 1000; i++) {
			sortObject.hash(partObject.locations);
			sortObject.kvSort();
			sortObject.findCellStartEnd();
			sortObject.reorder(partObject.locations, partObject.sortedLoc);
			partObject.countNeighbors(sortObject);
		}
		t = clock() - t;
		nnsTime = ((float)t) / CLOCKS_PER_SEC;
		printf("NNS time %0.3f\n", nnsTime);
	}

	// All-to-all
	{
		t = clock();
		for (int i = 0; i < 1000; i++) {
			partObject.countNeighborsN2(cellSize);
		}
		t = clock() - t;
		ataTime = ((float)t) / CLOCKS_PER_SEC;
		printf("All-to-all time %0.3f\n", ataTime);
	}

	printf("\nPerformance difference: %.1fx\n", ataTime / nnsTime);

#endif

	return 0;
}

/*

<Performance test output on desktop>

// NNS algorithm ran faster using the physical core count instead of logical core count
// Used an Intel i9-7920X with a base frequency of 2.9GHz

omp_get_max_threads() = 24
omp_get_max_threads() = 12

Running 1000 iterations of NNS and all-to-all

Simulation space is about: 60 x 60 x 60
Particle count:            3600

Average particles per non buffer cell: 2.08 (2.0 is not sparse)
At greater densities the NNS performance gain will decrease

NNS time 1.067
All-to-all time 34.745

Performance difference: 32.6x

*/