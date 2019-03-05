
find_package(VTK REQUIRED)
include ("${VTK_USE_FILE}")

ADD_EXECUTABLE( DummyVolumeGenerator DummyVolumeGenerator.cxx )
TARGET_LINK_LIBRARIES(DummyVolumeGenerator ${VTK_LIBRARIES} )

ADD_EXECUTABLE( VolumeTransform VolumeTransform.cxx )
TARGET_LINK_LIBRARIES(VolumeTransform ${VTK_LIBRARIES} )

ADD_EXECUTABLE( AverageVolumes AverageVolumes.cxx )
TARGET_LINK_LIBRARIES(AverageVolumes ${VTK_LIBRARIES} )
