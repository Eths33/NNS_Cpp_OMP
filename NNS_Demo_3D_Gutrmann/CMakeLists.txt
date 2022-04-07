cmake_minimum_required (VERSION 3.9)
project (nearest_neighbor_3D_search LANGUAGES CXX)

# Build type 
if( CMAKE_BUILD_TYPE STREQUAL "" )
	set( CMAKE_BUILD_TYPE "Debug" )
endif()

if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
	add_definitions(-D_DEBUG)
	
	message(" <<< Building Debug >>>")
    if (MSVC)
        # warning level 4
		# https://docs.microsoft.com/en-us/cpp/build/reference/compiler-option-warning-level?view=msvc-170
        add_compile_options(/DEBUG /Od /Zi /W4) #/WX all warnings as errors
    else()
        # lots of warnings
		# https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html
        add_compile_options(-g -O0 -Wall -Wextra -pedantic) #-Werror all warnings as errors
    endif()

else()

	message("<<< Building Relsease >>>")
    if (MSVC)
        add_compile_options(/O2) # (Maximize Speed)
    else()
        add_compile_options(-O3) # https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
    endif()

endif()

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