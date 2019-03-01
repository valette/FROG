cmake_minimum_required(VERSION 2.8)
project(FROG)

SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")

OPTION(BUILD_SHARED_LIBS "Build vtkOpenSURF with shared libraries." ON)

find_package(OpenCV REQUIRED)
find_package(VTK REQUIRED NO_MODULE)
include(${VTK_USE_FILE})

ADD_SUBDIRECTORY(match)
ADD_SUBDIRECTORY(registration)
