cmake_minimum_required(VERSION 2.8.3) 

project(DCoL)

message("CMAKE_SYSTEM_PROCESSOR: ${CMAKE_SYSTEM_PROCESSOR}")

# CHECK
if( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
    message( FATAL_ERROR "Please select another Build Directory" )
endif()
if( CMAKE_SOURCE_DIR MATCHES " " )
	message( "Your Source Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()
if( CMAKE_BINARY_DIR MATCHES " " )
	message( "Your Build Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()

#include(FindPkgConfig)

# Set output paths
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

#OpenMP
option (USE_OpenMP "Use OpenMP" OFF)
if(USE_OpenMP)
  find_package(OpenMP REQUIRED)
  if(OPENMP_FOUND)
    add_definitions(-DUSE_OpenMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

SET( srcs 
        DistanceFunctions/EuclideanFunction.cpp
		DistanceFunctions/NormalizedEuclideanFunction.cpp
		DistanceFunctions/StdWeightedEuclideanFunction.cpp
		DistanceFunctions/OverlapFunction.cpp
		DistanceFunctions/VDMFunction.cpp
		Date.cpp
		DateContainer.cpp
		StringTokenizer.cpp
		Utils.cpp
		Matrix.cpp
		Dataset.cpp
		ExtendedDataset.cpp
		ComplexityMeasures.cpp
		Statistics.cpp
		SMO.cpp
		InputOptions.cpp
		#Main.cpp
		#Main_mex_ori.cpp
		Main_mex.cpp
	)

add_executable( dcol ${srcs} )
target_link_libraries( dcol ${ALL_LIBS} )
