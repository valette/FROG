cmake_minimum_required(VERSION 3.30)
project(FROG)

SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")

OPTION(BUILD_SHARED_LIBS "Build FROG with shared libraries." ON)
OPTION(USE_UINT_FOR_POINT_IDS "Use unsigned ints for point ids" ON)

IF (USE_UINT_FOR_POINT_IDS)
	ADD_DEFINITIONS(-D INT_PTIDS)
ENDIF (USE_UINT_FOR_POINT_IDS)

find_package(OpenCV REQUIRED)
find_package(VTK REQUIRED COMPONENTS CommonSystem FiltersGeneral FiltersHybrid ImagingColor IOGeometry IOExport IOPLY )
include_directories(${VTK_INCLUDE_DIRS})

ADD_SUBDIRECTORY(match)
ADD_SUBDIRECTORY(registration)
ADD_SUBDIRECTORY(tools)
ADD_SUBDIRECTORY(vtkOpenSURF3D)
