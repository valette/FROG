
find_package(OpenCV REQUIRED)
find_package(VTK REQUIRED)

include ("${VTK_USE_FILE}")

find_package(OpenMP)

if (OPENMP_FOUND)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

ADD_EXECUTABLE( frog frog.cxx image.cxx imageGroup.cxx stats.cxx )
TARGET_LINK_LIBRARIES(frog  ${VTK_LIBRARIES} stdc++fs  )

