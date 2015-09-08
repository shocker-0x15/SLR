cmake_minimum_required(VERSION 3.3)
set(CMAKE_CONFIGURATION_TYPES Debug Release)
project(SLR)

option(USE_LIBCPP "Use libc++ instead of libstdc++." ON)

# macro (set_xcode_property TARGET XCODE_PROPERTY XCODE_VALUE)
# set_property (TARGET ${TARGET} PROPERTY XCODE_ATTRIBUTE_${XCODE_PROPERTY}
#     ${XCODE_VALUE})
# endmacro (set_xcode_property)

if(NOT MSVC)
	# C++11と標準ライブラリのサポートチェック
	include(CheckCXXCompilerFlag)
	CHECK_CXX_COMPILER_FLAG("-std=gnu++11" GNUPP11_SUPPORTED)
	CHECK_CXX_COMPILER_FLAG("-std=c++11" CPP11_SUPPORTED)
	CHECK_CXX_COMPILER_FLAG("-stdlib=libc++" LIBCPP_SUPPORTED)
	CHECK_CXX_COMPILER_FLAG("-stdlib=libstdc++" LIBSTDCPP_SUPPORTED)

	# GNU C++11を優先する。
	if(GNUPP11_SUPPORTED)
		set(CMAKE_CXX_FLAGS "-std=gnu++11")
	elseif(CPP11_SUPPORTED)
		set(CMAKE_CXX_FLAGS "-std=c++11")
	endif()

	# 標準ライブラリの選択	
	if(LIBCPP_SUPPORTED AND USE_LIBCPP)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
	else()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
	endif()
endif()

set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

set(EXTLIBS_Assimp_include "/usr/local/include/" CACHE PATH "Assimp include directory")
set(EXTLIBS_Assimp_lib "/usr/local/lib/" CACHE PATH "Assimp library directory")
set(EXTLIBS_libpng16_include "/usr/local/include/" CACHE PATH "libpng 1.6 include directory")
set(EXTLIBS_libpng16_lib "/usr/local/lib/" CACHE PATH "libpng library directory")
set(EXTLIBS_OpenEXR22_include "/usr/local/include/OpenEXR" CACHE PATH "OpenEXR 2.2 include directory")
set(EXTLIBS_OpenEXR22_lib "/usr/local/lib/" CACHE PATH "OpenEXR library directory")

set(include_dirs "${EXTLIBS_Assimp_include};${EXTLIBS_libpng16_include};${EXTLIBS_OpenEXR22_include}")
set(lib_dirs "${EXTLIBS_Assimp_lib};${EXTLIBS_libpng16_lib};${EXTLIBS_OpenEXR22_lib}")
set(linker_flags "-lAssimp -lpng -lHalf -lIlmImf-2_2")

if(msvc)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

file(GLOB SLR_SOURCES
	SLR/*.h
	SLR/*.cpp
	SLR/Core/*.h
	SLR/Core/*.cpp
	SLR/Accelerator/*.h
	SLR/Accelerator/*.cpp
	SLR/BasicTypes/*.h
	SLR/BasicTypes/*.cpp
	SLR/BSDFs/*.h
	SLR/BSDFs/*.cpp
	SLR/EDFs/*.h
	SLR/EDFs/*.cpp
	SLR/Cameras/*.h
	SLR/Cameras/*.cpp
	SLR/Helper/*.h
	SLR/Helper/*.cpp
	SLR/Memory/*.h
	SLR/Memory/*.c
	SLR/Memory/*.cpp
	SLR/Renderers/*.h
	SLR/Renderers/*.cpp
	SLR/SceneGraph/*.h
	SLR/SceneGraph/*.cpp
	SLR/Surface/*.h
	SLR/Surface/*.cpp
	SLR/SurfaceMaterials/*.h
	SLR/SurfaceMaterials/*.cpp
	SLR/Textures/*.h
	SLR/Textures/*.cpp
	)

source_group("Core" REGULAR_EXPRESSION "SLR/Core/.*\.(h|c|hpp|cpp)")
source_group("Accelerator" REGULAR_EXPRESSION "SLR/Accelerator/.*\.(h|c|hpp|cpp)")
source_group("BasicTypes" REGULAR_EXPRESSION "SLR/BasicTypes/.*\.(h|c|hpp|cpp)")
source_group("BSDFs" REGULAR_EXPRESSION "SLR/BSDFs/.*\.(h|c|hpp|cpp)")
source_group("EDFs" REGULAR_EXPRESSION "SLR/EDFs/.*\.(h|c|hpp|cpp)")
source_group("Cameras" REGULAR_EXPRESSION "SLR/Cameras/.*\.(h|c|hpp|cpp)")
source_group("Helper" REGULAR_EXPRESSION "SLR/Helper/.*\.(h|c|hpp|cpp)")
source_group("Memory" REGULAR_EXPRESSION "SLR/Memory/.*\.(h|c|hpp|cpp)")
source_group("Renderers" REGULAR_EXPRESSION "SLR/Renderers/.*\.(h|c|hpp|cpp)")
source_group("SceneGraph" REGULAR_EXPRESSION "SLR/SceneGraph/.*\.(h|c|hpp|cpp)")
source_group("Surface" REGULAR_EXPRESSION "SLR/Surface/.*\.(h|c|hpp|cpp)")
source_group("SurfaceMaterials" REGULAR_EXPRESSION "SLR/SurfaceMaterials/.*\.(h|c|hpp|cpp)")
source_group("Textures" REGULAR_EXPRESSION "SLR/Textures/.*\.(h|c|hpp|cpp)")

include_directories(${include_dirs})
# link_directories(${lib_dirs})
foreach(lib_dir ${lib_dirs})
	link_directories(${lib_dir})
endforeach()
add_executable(SLR ${SLR_SOURCES})
target_link_libraries(SLR ${linker_flags})
