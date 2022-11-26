
# find_package(VTK REQUIRED)
# include ("${VTK_USE_FILE}")

find_package( OpenMP REQUIRED )
set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

SET(TOOL_FILES
AverageVolumes
CheckDiffeomorphism
DummyVolumeGenerator
MeshTransform
PointsTransform
VolumeTransform
)

FOREACH(loop_var ${TOOL_FILES})
  ADD_EXECUTABLE( ${loop_var} ${loop_var}.cxx )
  TARGET_LINK_LIBRARIES( ${loop_var} ${VTK_LIBRARIES} )
ENDFOREACH(loop_var)
