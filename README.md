# Coarse Triangulations

Given an input image, compute a coarse triangulation approximation of the image. Supported approximations: piecewise constant, linear, quadratic per triangle. An initial triangulation on the image is provided by the TRIM algorithm, and then a gradient descent mesh moving method is used to further align the triangular mesh to the image to minimize approximation error.

## Requirements

For speed purposes, an Nvidia GPU is needed with CUDA support. Other installations not included in this repo include MATLAB and `ffmpeg` (for Linux; this is not strictly necessary and is just used for output purposes).

## Compiling

On Linux, compiling can be done by
```
mkdir build
cd build
cmake ..
make
```

This will compile the main file as well as the tests, which can be found in the `bin` directory of the build.

## Running Triangulations

Put an image (for the sake of example, `image.jpg`) to be triangulated in the `images` directory. Then from the `build` directory, run ```./CoarseTriangulation image.jpg``` This will output a screenshot of the initial triangulation, the final coarse triangulation, a video of all intermediate steps, and graphs of the energy function and total step over time. All these can be found in the `outputs` directory. Note that running the program multiple times will overwrite these files, which may require root access.
