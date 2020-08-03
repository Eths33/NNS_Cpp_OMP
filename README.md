# NNS_Cpp_OMP
Nearest Neighbor Search in 3D Space

Author: Greg Gutmann

Designed as an introduction to one implementation of a nearest neighbor search. Coded in a C/C++ style.
It was created as a teaching tool for my work at Tokyo Institute of Technology. 

In addition to using 3D space instead of 2D space like Eths33/NNS_Cpp_SingleThread, it also has some multi-threading. 
OpenMP used in: NNS::hash, Particle::countNeighbors, and Particle::countNeighborsN2

Simplicity was preferred over performance in this code.

# Main steps are as follows

1. hash              // Function that takes particle locations and ouputs the cell ID that they are in

2. kvSort            // Function that uses C++'s std::sort, can be replaced with any key value sort

3. findCellStartEnd  // Function to find the start and stop particle indexies for each cell

4. reorder           // Function to reorder the primary data structures that are used

5. countNeighbors    // Function that uses the results of the NNS to do the users work

# Building

## Linux 

cd NNS_Demo_Gutrmann

mkdir build

cd build

cmake ..

make

## Windows 

Download CMake GUI

Source directory is NNS_Demo_Gutrmann

Build direcotory is NNS_Demo_Gutrmann/build

<it will ask you if it should make the build folder if it doesn't exist>

Click Configure

Click Generate

Go into build folder and click Visual Studio project file
