### Building the project 

From the main folder that contains the CMakeLists.txt run the following commands in the terminal.

Debug:
```
mkdir Debug
cd Debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```

Release:
```
mkdir Release
cd Release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```