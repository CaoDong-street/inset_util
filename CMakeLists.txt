cmake_minimum_required(VERSION 2.8.3)
project(inset_util)

# set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-deprecated-declarations")

set(BUILD_SHARED_LIBS ON)

get_filename_component(INSET_UTIL_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(iris_DIR "${INSET_UTIL_DIR}/../../../install/lib/cmake/iris")
set(mosek_DIR "${INSET_UTIL_DIR}/../../../install/lib/cmake/mosek")
find_package(iris REQUIRED)
find_package(mosek REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(casadi REQUIRED)
include_directories(
    include
    ${EIGEN3_INCLUDE_DIRS}
    ${IRIS_INCLUDE_DIRS}
    ${MOSEK_INCLUDE_DIRS})
link_libraries("${IRIS_LIBRARY_DIRS}/libiris_geometry.so.0")
link_directories(${IRIS_LIBRARY_DIRS})

add_library(inset SHARED include/inset_util/inset.cpp include/inset_util/inset.h)
set_target_properties(inset PROPERTIES SOVERSION 0)
target_link_libraries(inset mosek casadi)

install(FILES "${PROJECT_NAME}Config.cmake" "${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION "share/${PROJECT_NAME}/cmake")

install(
    TARGETS inset
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

install(DIRECTORY include/inset_basis include/inset_geometry include/inset_util DESTINATION include)
