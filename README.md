# WifiVisualizer
Wi-Fi signal strength visualizer based on ray tracing/ray marching methods.
### Libraries used:
EasyBMP, OpenGL Mathematics (glm), TinyObjLoader 

### Compiling:
OpenMP is required to compile the program. It is included by default in modern compilers (gcc-7 and last versions of clang).

g++ -O2 -fopenmp -o main src/Main.cpp src/Scene.cpp src/Tracer.cpp src/EasyBMP/EasyBMP.cpp src/TinyObjLoader/tiny_obj_loader.cc

### Usage:
CLI. Specify config file as 1st param.

### Config file structure:
config1.txt with an explanation to each line
```
  1280 720            //output image resolution
  Flat.obj          //path to .obj file containing model of room/house in which the distribution of wi-fi rays is needed to be visualized.
  14000 1000 100      //router coords inside the model
  100000              //initial signal strength, adjust according to room size and wifi router strength
  12500 1000 25000    //camera coords
  0 0 -1              //camera "forward" vector - defines camera direction
  0 1 0               //camera "up" vector - defines camera incline
  90 74               //camera FOV (horizontal and vertical respectfully)
  FALSE               //turns autocontrast filter on/off (experimental, recomennded to be set to FALSE)
  FALSE               //turns median filtering on/off. (experimental, recomennded to be set to FALSE)
  ```
### Output:
  Resulting image is placed in /img folder with a name corresponding to config file that was used.
  
P.S This program is purely just a PoC, software renderer and ray marching test
