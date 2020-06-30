# Coarse Triangulations

Given an input image, compute a coarse triangulation approximation of the image. Supported approximations: piecewise constant, linear, quadratic per triangle. An initial triangulation on the image is provided by the TRIM algorithm, and then a gradient descent mesh moving method is used to further align the triangular mesh to the image to minimize approximation error.

## Requirements

For speed purposes, an Nvidia GPU is needed with CUDA support. Other installations not included in this repo include MATLAB and `ffmpeg` (for Linux; this is not strictly necessary and is just used for output purposes).

## Compiling

On Linux, compiling can be done by
```
mkdir build
cd build
cmake ../
cmake --build .
```

If any modifications are made to files in `src` or to `main.cpp` after that, running `make` in the `build` directory will recompile the main file.

In order to compile the test folders on a Linux system, the prerequisite `libgtest-dev` is needed. Move to the `src` folder and run the following:
```
g++ -c *.cpp
mv *.o ../tests
g++ -Wall -g -pthread *.o FILENAME -lgtest_main -lgtest -lpthread -o OUTPUT
```
where FILENAME, OUTPUT are replaced by the desired names.

## Running Triangulations

Put an image (for the sake of example, `image.jpg`) to be triangulated in the `images` directory. Then from the `build` directory, run ```./CoarseTriangulation image.jpg DENSITY``` where the second argument `DENSITY` is a positive number between 0 and 1 indicating the density of the triangulation. This will output the initial TRIM triangulation, the final coarse triangulation, and a video of the intermediate steps all placed in the `outputs` directory.
