
find_package(VTK REQUIRED)
include ("${VTK_USE_FILE}")

find_package( OpenMP REQUIRED )
set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

ADD_EXECUTABLE( CheckDiffeomorphism CheckDiffeomorphism.cxx )
TARGET_LINK_LIBRARIES(CheckDiffeomorphism ${VTK_LIBRARIES} )

ADD_EXECUTABLE( DummyVolumeGenerator DummyVolumeGenerator.cxx )
TARGET_LINK_LIBRARIES(DummyVolumeGenerator ${VTK_LIBRARIES} )

ADD_EXECUTABLE( VolumeTransform VolumeTransform.cxx )
TARGET_LINK_LIBRARIES(VolumeTransform ${VTK_LIBRARIES} )

ADD_EXECUTABLE( AverageVolumes AverageVolumes.cxx )
TARGET_LINK_LIBRARIES(AverageVolumes ${VTK_LIBRARIES} )
