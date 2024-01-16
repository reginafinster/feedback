cmake_minimum_required(VERSION 2.6)
project(qssemanybody)

#Bring the headers into the project (full or relative path)
include_directories(itensor)
#add directory for own classes 
# include_directories(include)
#Set ITensor directory
set(ITENSOR_DIR /home/droenner/itensor/)

#Bring the headers into the project (full or relative path)
include_directories(${ITENSOR_DIR}
        ${ITENSOR_DIR}itensor/
        ${ITENSOR_DIR}itensor/detail
        ${ITENSOR_DIR}itensor/itdata
        ${ITENSOR_DIR}itensor/mps
        ${ITENSOR_DIR}itensor/mps/lattice/
        ${ITENSOR_DIR}itensor/mps/sites/
        ${ITENSOR_DIR}itensor/tensor
        ${ITENSOR_DIR}itensor/util
        )

#Link the Itensor library (NOTE: Compilation now works fine without these two lines)
# add_library(itensor STATIC IMPORTED)
# set_property(TARGET itensor PROPERTY IMPORTED_LOCATION ${ITENSOR_DIR}lib/libitensor-g.a)

#Set a variable with all the new flags
set(ITENSOR_FLAGS "-DPLATFORM_lapack -D__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES=0")
set(ITENSOR_LINK_FLAGS "-DPLATFORM_lapack -L/home/droenner/itensor/lib -litensor -lpthread -L/usr/lib -lblas -llapack")

#Append the new flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${ITENSOR_FLAGS}")
file(GLOB SOURCES "classes/*.cpp")

add_executable(qssemanybody main.cpp ${SOURCES})
target_link_libraries(qssemanybody itensor "${ITENSOR_LINK_FLAGS}")


