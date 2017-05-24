# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( optimal_convex_partition )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )
set(MAIN_FILE "optimal_convex_partition.cpp")

if ( CGAL_FOUND )
# CGAL static library creation
	include( ${CGAL_USE_FILE} )
	include( ${CMAKE_CURRENT_LIST_DIR}/CGAL_CreateSingleSourceCGALLib.cmake )
	set(CMAKE_CXX_STANDARD 11)
	create_single_source_cgal_lib( ${MAIN_FILE} )

	#Matlab and mex configurations
	set(MATLAB_ROOT "~/opt/MATLAB/R2016b")
	set(CMAKE_CXX_STANDARD 11)
	set(MEX_SAVE_PATH "${CMAKE_SOURCE_DIR}")
	set(FLAGS "-fPIC -DMX_COMPAT_32 -DMATLAB_MEX_FILE")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 ${FLAGS}")
	set(WRAP "mex_wrap")
	set(WRAP_FILE "${WRAP}.cpp")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")

	#include mex.h dir
	include_directories("${MATLAB_ROOT}/extern/include")
	# add_library(${PROJECT_NAME} ${MAIN_FILE})
	add_library(${WRAP} OBJECT ${WRAP_FILE})
	ADD_CUSTOM_COMMAND(TARGET ${PROJECT_NAME}
	    POST_BUILD
	    COMMAND mex  -L./  -l${PROJECT_NAME}  -output ${PROJECT_NAME} ${CMAKE_BINARY_DIR}/CMakeFiles/${WRAP}.dir/${WRAP_FILE}.o
)


else()
	message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

