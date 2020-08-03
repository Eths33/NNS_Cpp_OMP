cmake_minimum_required (VERSION 3.9)
project (nearest_neighbor_3D_search LANGUAGES CXX)

file( GLOB HEADER_FILES "header/*.hpp" )
file( GLOB SOURCE_FILES "source/*.cpp" )

source_group( "Header" FILES ${HEADER_FILES} )
source_group( "Source" FILES ${SOURCE_FILES} )

set( ALL_SAMPLE_FILES ${HEADER_FILES} ${SOURCE_FILES} )

add_executable (nearest_neighbor_3D_search ${ALL_SAMPLE_FILES})

target_compile_features(nearest_neighbor_3D_search PUBLIC cxx_std_17)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(nearest_neighbor_3D_search PUBLIC OpenMP::OpenMP_CXX)
endif()

IF (WIN32)
list(APPEND CMAKE_VS_SDK_INCLUDE_DIRECTORIES "$(VC_IncludePath);$(WindowsSDK_IncludePath)")
list(APPEND CMAKE_VS_SDK_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/header")
ELSE()
include_directories( "header" )
ENDIF()