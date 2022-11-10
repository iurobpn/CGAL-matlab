# CGAL matlab
## Warning
I have updated the CMakeLists.txt file to use the cmake FindMatlab module in the optimal_convex_partition_2 function. I will be updating the other function soon. This will make it compatible with CGAL 5.x and above.

## Descripion
CGAL functions in matlab using mex tool and cmake for compilation. 

three functions are provided at the moment:

```
optimal_convex_partition_2
greene_approx_convex_partition_2
convex_decomposition_3
```

The first two are for partition of a 2D nonconvex polygon into convex polygons. The last one is for decomposition of 3D polyhedron into convex polyhedrons. See example on how to use the functions in matlab. 

You can use these functions as a template to add another ones, please share if you adapt new cgal functions.


## Install
To install on ubuntu (tested on ubuntu 20.04) first install cgal(below version 5.0) and boost from the ubuntu repositories:

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

You have to provide the matlab path or edit the MATLAB_ROOT variable in the CMakeLists.txt file of the function you weant to compile.

Since version 5.0, CGAL is a header only library and this may change the compilation procedure. I did not tested with version 5.0 and above yet. This examples are tested with Matlab R2019b, cmake 3.17.3, and cgal 4.11 from Ubuntu 18.04 repositories (libcgal-dev package).

