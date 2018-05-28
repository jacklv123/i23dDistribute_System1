CMAKE_MINIMUM_REQUIRED(VERSION 2.6.2)
PROJECT(I23D)

# Add modules
# Find dependencies:
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/i23dMVS/build/Modules)

# Define helper functions and macros.
cmake_policy(SET CMP0011 OLD)
INCLUDE(i23dMVS/build/Utils.cmake)
if(ENABLE_PRECOMPILED_HEADERS)
	INCLUDE(i23dMVS/build/Cotire.cmake)
endif()

FIND_PACKAGE(Eigen ${SYSTEM_PACKAGE_REQUIRED})
if(EIGEN_FOUND)
	include_directories(${EIGEN_INCLUDE_DIRS})
	add_definitions(${EIGEN_DEFINITIONS} -D_USE_EIGEN)
endif()

ADD_SUBDIRECTORY(i23dMVS)
ADD_SUBDIRECTORY(i23dSFM)
ADD_SUBDIRECTORY(i23dDistribute_System)

