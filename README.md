# CGAL-matlab
CGAL functions in matlab using mex tool e cmake for compilation 

To build, run:

```{bash}
mkdir build
cd build
cmake -DCGAL_DIR=$CGAL_DIR
```

where '$CGAL_DIR' is the directory where you compiled cgal. You can also use the ubuntu package libcgal-dev and drop the '-DCGAL_DIR=$CGAL_DIR' after cmake. 

Since version 5.0, CGAL is a header only library and this may change the compilation procedure. I did not tested with version 5.0 and above yet. This examples are tested with Matlab R2019b, cmake 3.17.3, and cgal 4.11 from Ubuntu 18.04 repositories (libcgal-dev package).


