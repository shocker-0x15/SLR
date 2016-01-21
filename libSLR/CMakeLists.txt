set(include_dirs "${EXTLIBS_OpenEXR22_include}")
set(lib_dirs "${EXTLIBS_OpenEXR22_lib}")
set(linker_flags "-lHalf")

file(GLOB libSLR_Sources
	*.h
	*.cpp
	Core/*.h
	Core/*.cpp
	Accelerator/*.h
	Accelerator/*.cpp
	BasicTypes/*.h
	BasicTypes/*.cpp
	BSDFs/*.h
	BSDFs/*.cpp
	EDFs/*.h
	EDFs/*.cpp
	Cameras/*.h
	Cameras/*.cpp
	Helper/*.h
	Helper/*.cpp
	Memory/*.h
	Memory/*.c
	Memory/*.cpp
	Renderers/*.h
	Renderers/*.cpp
	SceneGraph/*.h
	SceneGraph/*.cpp
	Surface/*.h
	Surface/*.cpp
	SurfaceMaterials/*.h
	SurfaceMaterials/*.cpp
	Textures/*.h
	Textures/*.cpp
	)
file(GLOB libSLR_Sources_excluded
	 Memory/StackAllocator.*
     Memory/MSpaceAllocator.*
     Memory/dlmalloc.*
    )
list(REMOVE_ITEM libSLR_Sources ${libSLR_Sources_excluded})

source_group("SLR" REGULAR_EXPRESSION ".*\.(h|c|hpp|cpp)")
source_group("SLR\\Core" REGULAR_EXPRESSION "Core/.*\.(h|c|hpp|cpp)")
source_group("SLR\\Accelerator" REGULAR_EXPRESSION "Accelerator/.*\.(h|c|hpp|cpp)")
source_group("SLR\\BasicTypes" REGULAR_EXPRESSION "BasicTypes/.*\.(h|c|hpp|cpp)")
source_group("SLR\\BSDFs" REGULAR_EXPRESSION "BSDFs/.*\.(h|c|hpp|cpp)")
source_group("SLR\\EDFs" REGULAR_EXPRESSION "EDFs/.*\.(h|c|hpp|cpp)")
source_group("SLR\\Cameras" REGULAR_EXPRESSION "Cameras/.*\.(h|c|hpp|cpp)")
source_group("SLR\\Helper" REGULAR_EXPRESSION "Helper/.*\.(h|c|hpp|cpp)")
source_group("SLR\\Memory" REGULAR_EXPRESSION "Memory/.*\.(h|c|hpp|cpp)")
source_group("SLR\\Renderers" REGULAR_EXPRESSION "Renderers/.*\.(h|c|hpp|cpp)")
source_group("SLR\\Surface" REGULAR_EXPRESSION "Surface/.*\.(h|c|hpp|cpp)")
source_group("SLR\\SurfaceMaterials" REGULAR_EXPRESSION "SurfaceMaterials/.*\.(h|c|hpp|cpp)")
source_group("SLR\\Textures" REGULAR_EXPRESSION "Textures/.*\.(h|c|hpp|cpp)")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

include_directories(${include_dirs})
# link_directories(${lib_dirs})
foreach(lib_dir ${lib_dirs})
	link_directories(${lib_dir})
endforeach()
add_library(SLR SHARED ${libSLR_Sources})
target_link_libraries(SLR ${linker_flags})

install(TARGETS SLR CONFIGURATIONS Debug DESTINATION "${CMAKE_BINARY_DIR}/bin/Debug")
install(TARGETS SLR CONFIGURATIONS Release DESTINATION "${CMAKE_BINARY_DIR}/bin/Release")
