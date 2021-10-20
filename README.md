[#](#) CGAL-matlab
CGAL functions in matlab using mex tool e cmake for compilation. To install on ubuntu (tested on ubuntu 20.04) first install cgal(below version 5.0) and boost from the ubuntu repositories:

```{bash}
sudo apt install libboost-dev libcgal-dev
```

To compile one of the functions, say optimal convex partition 2:

```{bash}
cd path-to-repo/optimal_convex_partition_2/
mkdir build
cd build
cmake .. -DMATLAB_ROOT=path-to-matlab
make
```

You have to provide the matlab path or edit the MATLAB_ROOT file in the CMakeLists.txt file of the function you weant to compile.

Since version 5.0, CGAL is a header only library and this may change the compilation procedure. I did not tested with version 5.0 and above yet. This examples are tested with Matlab R2019b, cmake 3.17.3, and cgal 4.11 from Ubuntu 18.04 repositories (libcgal-dev package).

