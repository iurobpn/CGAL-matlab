# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.24)
project( convex_decomposition_3 )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(Matlab_ROOT_DIR $ENV{HOME}/opt/MATLAB/R2022b/)

find_package(Matlab REQUIRED)
# find_package(Armadillo REQUIRED)
find_package(CGAL REQUIRED COMPONENTS Core)

#FindMatlab()
include_directories(/usr/local/include .)
matlab_add_mex(
	NAME convex_decomposition_3
	SRC convex_decomposition_3.cpp
	LINK_TO CGAL::CGAL CGAL::CGAL_Core)

