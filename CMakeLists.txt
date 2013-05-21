cmake_minimum_required(VERSION 2.8)
project(Estimation CXX)

find_package(Boost REQUIRED)

set(ESTIMATION_SRC_DIRS ${PROJECT_SOURCE_DIR}/src)
set(ESTIMATION_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)
# # currently only building certain files.
# set(ESTIMATION_SRCS
# )
# set(ESTIMATION_HEADERS
#     common.h
#     invert_matrix.hpp
#     # expm.hpp
# )

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(${ESTIMATION_INCLUDE_DIRS})

# Boost
link_directories(${Boost_LIBRARY_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

add_library(Estimation ${ESTIMATION_INCLUDE_DIRS} ${ESTIMATION_SRC_DIRS})
set_target_properties(Estimation PROPERTIES LINKER_LANGUAGE CXX)

install(
    TARGETS Estimation
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

############################################

# PROJECT(gCommon)

# SET(GSTATESPACELIB_SRCS
# gKalmanFilter.h
# gKalmanFilter.cpp
# gControl.h
# gControl.cpp
# gSSCommon.h
# gProbability.h
# gProbability.cpp
# InvertMatrix.h
# UKF.cpp
# UKF.h
# matrix2.h
# EKF.cpp
# EKF.h
# )


# SET( GSTATESPACELIB_HEADERS
# gKalmanFilter.h
# gControl.h
# gSSCommon.h
# InvertMatrix.h
# gProbability.h
# UKF.h
# EKF.h
# matrix2.h
# )

# ADD_LIBRARY(gStateSpace ${GSTATESPACELIB_SRCS})

# # install headers
# INSTALL( FILES ${GSTATESPACELIB_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/gStateSpaceLib )

# # install libraries
# INSTALL(TARGETS
# gStateSpace
# RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}
# LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}
# ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}
# )
