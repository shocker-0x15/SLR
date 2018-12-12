cmake_minimum_required(VERSION 3.3)
set(CMAKE_CONFIGURATION_TYPES Debug Release)
project(SLR)

option(USE_LIBCPP "Use libc++ instead of libstdc++." ON)

# macro (set_xcode_property TARGET XCODE_PROPERTY XCODE_VALUE)
# set_property (TARGET ${TARGET} PROPERTY XCODE_ATTRIBUTE_${XCODE_PROPERTY}
#     ${XCODE_VALUE})
# endmacro (set_xcode_property)

if(MSVC)
    add_definitions(/MP)

    set(CMAKE_C_FLAGS_DEBUG "/Od /DDEBUG /ZI")
    set(CMAKE_CXX_FLAGS_DEBUG "/Od /DDEBUG /ZI")
    set(CMAKE_C_FLAGS_RELEASE "/O2 /DNDEBUG /Zi")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG /Zi")

    set_property(GLOBAL APPEND PROPERTY LINK_FLAGS_DEBUG /DEBUG:FASTLINK)
    set_property(GLOBAL APPEND PROPERTY LINK_FLAGS_RELEASE /DEBUG:FASTLINK)
else()
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

    set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif()

set(EXTLIBS_Assimp_include "/usr/local/include/" CACHE PATH "Assimp include directory")
set(EXTLIBS_Assimp_lib "/usr/local/lib/" CACHE PATH "Assimp library directory")
set(EXTLIBS_libpng16_include "/usr/local/include/" CACHE PATH "libpng 1.6 include directory")
set(EXTLIBS_libpng16_lib "/usr/local/lib/" CACHE PATH "libpng library directory")
set(EXTLIBS_OpenEXR22_include "/usr/local/include/OpenEXR" CACHE PATH "OpenEXR 2.2 include directory")
set(EXTLIBS_OpenEXR22_lib "/usr/local/lib/" CACHE PATH "OpenEXR library directory")

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_definitions(-DOPENEXR_DLL)
endif()

# OS Xにおけるrun path処理の有効化
set(CMAKE_MACOSX_RPATH 1)

# 各プロジェクトのCMakeLists.txtを呼び出す。
add_subdirectory(libSLR)
add_subdirectory(libSLRSceneGraph)
add_subdirectory(HostProgram)

# ビルド依存関係を設定
add_dependencies(SLRSceneGraph SLR)
add_dependencies(HostProgram SLR SLRSceneGraph)
