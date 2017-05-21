set(include_dirs "${EXTLIBS_Assimp_include};${EXTLIBS_libpng16_include};${EXTLIBS_OpenEXR22_include};${CMAKE_SOURCE_DIR}")
set(lib_dirs "${EXTLIBS_Assimp_lib};${EXTLIBS_libpng16_lib};${EXTLIBS_OpenEXR22_lib};")
set(libs "Assimp;png;Half;IlmImf;SLR")

file(GLOB libSLRSceneGraph_Sources
     *.h
     *.cpp
     Scene/*.h
     Scene/*.cpp
     Parser/BuiltinFunctions/*.h
     Parser/BuiltinFunctions/*.cpp
     Parser/*.h
     Parser/*.cpp
     Parser/*.hh
     Parser/*.cc
     Helper/*.h
     Helper/*.cpp
    )
file(GLOB libSLRSceneGraph_Sources_excluded
     Parser/main.cpp
    )
# workaround for CMake 3.4.1 on Windows.
set(libSLRSceneGraph_Sources_excluded "${libSLRSceneGraph_Sources_excluded} dummy")
list(REMOVE_ITEM libSLRSceneGraph_Sources ${libSLRSceneGraph_Sources_excluded})

source_group(""         REGULAR_EXPRESSION        ".*\.(h|c|hpp|cpp)")
source_group("Scene"  REGULAR_EXPRESSION  "Scene/.*\.(h|c|hpp|cpp)")
source_group("Parser" REGULAR_EXPRESSION "Parser/.*\.(h|c|hh|cc|hpp|cpp)")
source_group("Parser\\Builtin Functions" REGULAR_EXPRESSION "Parser/BuiltinFunctions/.*\.(h|c|hpp|cpp)")
source_group("Helper" REGULAR_EXPRESSION "Helper/.*\.(h|c|hpp|cpp)")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(MSVC)
    add_definitions(-DSLR_SCENEGRAPH_API_EXPORTS)
endif()

include_directories(${include_dirs})
# link_directories(${lib_dirs})
foreach(lib_dir ${lib_dirs})
    link_directories(${lib_dir})
endforeach()
add_library(SLRSceneGraph SHARED ${libSLRSceneGraph_Sources})
foreach(lib ${libs})
    target_link_libraries(SLRSceneGraph PRIVATE ${lib})
endforeach()

set_target_properties(SLRSceneGraph PROPERTIES INSTALL_RPATH "@executable_path")
install(TARGETS SLRSceneGraph CONFIGURATIONS Debug DESTINATION "${CMAKE_BINARY_DIR}/bin/Debug")
install(TARGETS SLRSceneGraph CONFIGURATIONS Release DESTINATION "${CMAKE_BINARY_DIR}/bin/Release")
