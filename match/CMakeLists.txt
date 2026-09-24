find_package( OpenMP REQUIRED )
find_package( Boost COMPONENTS filesystem REQUIRED  iostreams REQUIRED )

option(USE_SSE_FOR_MATCHING "Use SSE/AVX instructions for matching (may crash on old computers)" OFF)

set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

if(USE_SSE_FOR_MATCHING)
  add_definitions(-D USE_SSE_FOR_MATCHING)
  if(MSVC)
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX" )
  else()
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx" )
  endif()
endif(USE_SSE_FOR_MATCHING)

include_directories(${Boost_INCLUDE_DIRS})
ADD_EXECUTABLE( match match.cpp )
target_link_libraries(match ${Boost_LIBRARIES} ${VTK_LIBRARIES})
