set(include_dirs "${EXTLIBS_OpenEXR22_include};${CMAKE_SOURCE_DIR}")
set(lib_dirs "")
set(libs "SLR;SLRSceneGraph")

file(GLOB HostProgram_Sources
     main.cpp
     StopWatch.h
    )

source_group("" REGULAR_EXPRESSION ".*\.(h|c|hpp|cpp)")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

include_directories(${include_dirs})
# link_directories(${lib_dirs})
foreach(lib_dir ${lib_dirs})
    link_directories(${lib_dir})
endforeach()
add_executable(HostProgram ${HostProgram_Sources})
foreach(lib ${libs})
    target_link_libraries(HostProgram PRIVATE ${lib})
endforeach()

set_target_properties(HostProgram PROPERTIES INSTALL_RPATH "@executable_path")
install(TARGETS HostProgram CONFIGURATIONS Debug DESTINATION "${CMAKE_BINARY_DIR}/bin/Debug")
install(TARGETS HostProgram CONFIGURATIONS Release DESTINATION "${CMAKE_BINARY_DIR}/bin/Release")
