# Tests CMakeLists.txt

# test_stats
add_executable(test_stats test_stats.cxx ${CMAKE_SOURCE_DIR}/registration/stats.cxx)
target_include_directories(test_stats PRIVATE ${CMAKE_SOURCE_DIR}/registration)
add_test(NAME test_stats COMMAND test_stats)

# test_transform_io
add_executable(test_transform_io test_transform_io.cxx)
target_include_directories(test_transform_io PRIVATE
	${CMAKE_SOURCE_DIR}/tools
	${CMAKE_SOURCE_DIR}/vtkOpenSURF3D
	${VTK_INCLUDE_DIRS}
)
target_link_libraries(test_transform_io PRIVATE ${VTK_LIBRARIES})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
	target_link_libraries(test_transform_io PRIVATE stdc++fs)
endif()
add_test(NAME test_transform_io COMMAND test_transform_io)
