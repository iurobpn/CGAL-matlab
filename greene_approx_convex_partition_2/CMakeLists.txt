# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( greene_approx_convex_partition)

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )
set(MAIN_FILE "greene_approx_convex_partition.cpp")
set(CMAKE_BUILD_TYPE "Release")
# change this variable to set matlab path
# set(MATLAB_ROOT "/usr/local/MATLAB/R2019b")

if ( CGAL_FOUND )
# CGAL static library creation
	include( ${CGAL_USE_FILE} )
	include( ${CMAKE_CURRENT_LIST_DIR}/CGAL_CreateSingleSourceCGALLib.cmake )
	set(CMAKE_CXX_STANDARD 11)
	create_single_source_cgal_lib( ${MAIN_FILE} )

	#Matlab and mex configurations
	set(CMAKE_CXX_STANDARD 11)
	set(MEX_SAVE_PATH "${CMAKE_SOURCE_DIR}")
	set(FLAGS "-fPIC -DMX_COMPAT_32 -DMATLAB_MEX_FILE")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 ${FLAGS}")
	set(WRAP "mex_wrap")
	set(WRAP_FILE "${WRAP}.cpp")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
	set(CGAL_LINK_FLAGS "-lboost_system -lpthread -lboost_thread  -lgmp -lCGAL -lCGAL_Core")

	#include mex.h dir
	include_directories("${MATLAB_ROOT}/extern/include")
	link_directories("/usr/lib/x86_64-linux-gnu")
	# add_library(${PROJECT_NAME} ${MAIN_FILE})
	add_library(${WRAP} OBJECT ${WRAP_FILE})
	ADD_CUSTOM_COMMAND(TARGET ${PROJECT_NAME}
	    POST_BUILD
	    COMMAND ${MATLAB_ROOT}/bin/mex  -L./  -outdir .. -output ${PROJECT_NAME} $<TARGET_OBJECTS:${WRAP}>  -l${PROJECT_NAME} -lboost_system -lpthread -lboost_thread  -lgmp -lCGAL -lCGAL_Core
)

	    # COMMAND mex  -L./  -output ${PROJECT_NAME} ${CMAKE_BINARY_DIR}/CMakeFiles/${WRAP}.dir/${WRAP_FILE}.o  -l${PROJECT_NAME} ${CGAL_LINK_FLAGS}

else()
	message(STATUS "This program requires the CGAL library, and will not be compiled.")
  # -lCGAL_Core
endif()

# message("CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS: ${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS}") 
# message("CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG: ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}") 
# message("CGAL_LIBRARIES: ${CGAL_LIBRARIES}") 
# message("CGAL_3RD_PARTY_LIBRARIES: ${CGAL_3RD_PARTY_LIBRARIES}") 
# message("CGAL_EXECUTABLE_TARGETS: ${CGAL_EXECUTABLE_TARGETS}") 
# message("CGAL_LIBRARY_TARGETS: ${CGAL_LIBRARY_TARGETS}") 
