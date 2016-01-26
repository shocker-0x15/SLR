set(include_dirs "${EXTLIBS_Assimp_include};${EXTLIBS_libpng16_include};${EXTLIBS_OpenEXR22_include};${CMAKE_SOURCE_DIR}")
set(lib_dirs "${EXTLIBS_Assimp_lib};${EXTLIBS_libpng16_lib};${EXTLIBS_OpenEXR22_lib};")
set(linker_flags "Assimp png Half IlmImf SLR")

file(GLOB libSLRSceneGraph_Sources
	*.h
	*.hpp
	*.cpp
	Helper/*.h
	Helper/*.cpp
	Parser/*.h
	Parser/*.hh
	Parser/*.cpp
	Parser/*.cc
	)
file(GLOB libSLRSceneGraph_Sources_excluded
	 Parser/main.cpp
    )
# workaround for CMake 3.4.1 on Windows.
set(libSLRSceneGraph_Sources_excluded "${libSLRSceneGraph_Sources_excluded} dummy")
list(REMOVE_ITEM libSLRSceneGraph_Sources ${libSLRSceneGraph_Sources_excluded})

source_group("SLR Scene Graph" REGULAR_EXPRESSION ".*\.(h|c|hpp|cpp)")
source_group("SLR Scene Graph\\Helper" REGULAR_EXPRESSION "Helper/.*\.(h|c|hpp|cpp)")
source_group("SLR Scene Graph\\Parser" REGULAR_EXPRESSION "Parser/.*\.(h|c|hh|cc|hpp|cpp)")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(MSVC)
	add_definitions(-DSLR_SCENEGRAPH_API_EXPORTS)
endif()

include_directories(${include_dirs})
# link_directories(${lib_dirs})
foreach(lib_dir ${lib_dirs})
	link_directories(${lib_dir})
endforeach()
add_library(SLRSceneGraph SHARED ${libSLRSceneGraph_Sources})
target_link_libraries(SLRSceneGraph ${linker_flags})

set_target_properties(SLRSceneGraph PROPERTIES INSTALL_RPATH "@executable_path")
install(TARGETS SLRSceneGraph CONFIGURATIONS Debug DESTINATION "${CMAKE_BINARY_DIR}/bin/Debug")
install(TARGETS SLRSceneGraph CONFIGURATIONS Release DESTINATION "${CMAKE_BINARY_DIR}/bin/Release")
