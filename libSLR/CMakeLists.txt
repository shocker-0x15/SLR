set(include_dirs "${EXTLIBS_OpenEXR22_include}")
set(lib_dirs "${EXTLIBS_OpenEXR22_lib}")
set(libs "Half")

file(GLOB libSLR_Sources
     *.h
     *.cpp
     External/AnalyticSkyDomeRadiance/*.h
     External/AnalyticSkyDomeRadiance/*.c
     MemoryAllocators/*.h
     MemoryAllocators/*.cpp
     BasicTypes/*.h
     BasicTypes/*.cpp
     Core/*.h
     Core/*.cpp
     RNG/*.h
     RNG/*.cpp
     SurfaceShape/*.h
     SurfaceShape/*.cpp
     MediumDistribution/*.h
     MediumDistribution/*.cpp
     BSDF/*.h
     BSDF/*.cpp
     EDF/*.h
     EDF/*.cpp
     PhaseFunction/*.h
     PhaseFunction/*.cpp
     Texture/*.h
     Texture/*.cpp
     SurfaceMaterial/*.h
     SurfaceMaterial/*.cpp
     MediumMaterial/*.h
     MediumMaterial/*.cpp
     Accelerator/*.h
     Accelerator/*.cpp
     Camera/*.h
     Camera/*.cpp
     Scene/*.h
     Scene/*.cpp
     Renderer/*.h
     Renderer/*.cpp
     Helper/*.h
     Helper/*.cpp
    )
file(GLOB libSLR_Sources_excluded
     MemoryAllocators/StackAllocator.*
     MemoryAllocators/MSpaceAllocator.*
     MemoryAllocators/dlmalloc.*
     Renderer/AMCMCPPMRenderer.*
    )
list(REMOVE_ITEM libSLR_Sources ${libSLR_Sources_excluded})

source_group(""                    REGULAR_EXPRESSION                    ".*\.(h|c|hpp|cpp)")
source_group("External\\Analytic Sky Dome Radiance"   REGULAR_EXPRESSION   "External/AnalyticSkyDomeRadiance/.*\.(h|c|hpp|cpp)")
source_group("Memory Allocators"   REGULAR_EXPRESSION   "MemoryAllocators/.*\.(h|c|hpp|cpp)")
source_group("Basic Types"         REGULAR_EXPRESSION         "BasicTypes/.*\.(h|c|hpp|cpp)")
source_group("Core"                REGULAR_EXPRESSION               "Core/.*\.(h|c|hpp|cpp)")
source_group("RNG"                 REGULAR_EXPRESSION                "RNG/.*\.(h|c|hpp|cpp)")
source_group("Surface Shape"       REGULAR_EXPRESSION       "SurfaceShape/.*\.(h|c|hpp|cpp)")
source_group("Medium Distribution" REGULAR_EXPRESSION "MediumDistribution/.*\.(h|c|hpp|cpp)")
source_group("BSDF"                REGULAR_EXPRESSION               "BSDF/.*\.(h|c|hpp|cpp)")
source_group("EDF"                 REGULAR_EXPRESSION                "EDF/.*\.(h|c|hpp|cpp)")
source_group("Phase Function"      REGULAR_EXPRESSION      "PhaseFunction/.*\.(h|c|hpp|cpp)")
source_group("Texture"             REGULAR_EXPRESSION            "Texture/.*\.(h|c|hpp|cpp)")
source_group("Surface Material"    REGULAR_EXPRESSION    "SurfaceMaterial/.*\.(h|c|hpp|cpp)")
source_group("Medium Material"     REGULAR_EXPRESSION     "MediumMaterial/.*\.(h|c|hpp|cpp)")
source_group("Accelerator"         REGULAR_EXPRESSION        "Accelerator/.*\.(h|c|hpp|cpp)")
source_group("Camera"              REGULAR_EXPRESSION             "Camera/.*\.(h|c|hpp|cpp)")
source_group("Scene"               REGULAR_EXPRESSION              "Scene/.*\.(h|c|hpp|cpp)")
source_group("Renderer"            REGULAR_EXPRESSION           "Renderer/.*\.(h|c|hpp|cpp)")
source_group("Helper"              REGULAR_EXPRESSION             "Helper/.*\.(h|c|hpp|cpp)")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(MSVC)
    add_definitions(-DSLR_API_EXPORTS)
endif()

include_directories(${include_dirs})
# link_directories(${lib_dirs})
foreach(lib_dir ${lib_dirs})
    link_directories(${lib_dir})
endforeach()
add_library(SLR SHARED ${libSLR_Sources})
foreach(lib ${libs})
    target_link_libraries(SLR PRIVATE ${lib})
endforeach()

install(TARGETS SLR CONFIGURATIONS Debug DESTINATION "${CMAKE_BINARY_DIR}/bin/Debug")
install(TARGETS SLR CONFIGURATIONS Release DESTINATION "${CMAKE_BINARY_DIR}/bin/Release")
