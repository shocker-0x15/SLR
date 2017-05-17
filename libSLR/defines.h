//
//  defines.h
//
//  Created by 渡部 心 on 2015/07/08.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_defines__
#define __SLR_defines__

// Platform defines
#if defined(_WIN32) || defined(_WIN64)
#   define SLR_Platform_Windows
#   if defined(__MINGW32__) // Defined for both 32 bit/64 bit MinGW
#       define SLR_Platform_Windows_MinGW
#   elif defined(_MSC_VER)
#       define SLR_Platform_Windows_MSVC
#   endif
#elif defined(__linux__)
#   define SLR_Platform_Linux
#elif defined(__APPLE__)
#   define SLR_Platform_OS_X
#elif defined(__OpenBSD__)
#   define SLR_Platform_OpenBSD
#endif

#ifdef SLR_Platform_Windows_MSVC
#   define NOMINMAX
#   define _USE_MATH_DEFINES
#   ifdef SLR_API_EXPORTS
#       define SLR_API __declspec(dllexport)
#   else
#       define SLR_API __declspec(dllimport)
#   endif
// MSVC 19.0 (Visual Studio 2015 Update 1) seems to have a problem related to a constexpr constructor.
#   define CONSTEXPR_CONSTRUCTOR
#   include <Windows.h>
#   undef near
#   undef far
#else
#   define SLR_API
#   define CONSTEXPR_CONSTRUCTOR constexpr
#endif

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstdarg>
#include <cmath>
#include <cfloat>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include <array>
#include <vector>
#include <deque>
#include <map>
#include <set>
#include <stack>

#include <chrono>
#include <limits>
#include <algorithm>
#include <memory>
#include <functional>

#ifdef DEBUG
#   define ENABLE_ASSERT
#endif

#ifdef SLR_Platform_Windows_MSVC
SLR_API void debugPrintf(const char* fmt, ...);
#else
#   define debugPrintf(fmt, ...) printf(fmt, ##__VA_ARGS__);
#endif

#ifdef ENABLE_ASSERT
#   define SLRAssert(expr, fmt, ...) if (!(expr)) { debugPrintf("%s @%s: %u:\n", #expr, __FILE__, __LINE__); debugPrintf(fmt"\n", ##__VA_ARGS__); abort(); } 0
#else
#   define SLRAssert(expr, fmt, ...)
#endif

#define SLRAssert_ShouldNotBeCalled() SLRAssert(false, "Should not be called!")
#define SLRAssert_NotImplemented() SLRAssert(false, "Not implemented yet!")

#define SLR_Minimum_Machine_Alignment 16
#define SLR_L1_Cacheline_Size 64

// For memalign, free, alignof
#if defined(SLR_Platform_Windows_MSVC)
#   include <malloc.h>
#   define SLR_memalign(size, alignment) _aligned_malloc(size, alignment)
#   define SLR_freealign(ptr) _aligned_free(ptr)
#   define SLR_alignof(T) __alignof(T)
#elif defined(SLR_Platform_OS_X) || defined(SLR_Platform_OpenBSD)
inline void* SLR_memalign(size_t size, size_t alignment) {
    void* ptr;
    if (posix_memalign(&ptr, alignment, size))
        ptr = nullptr;
    return ptr;
}
#   define SLR_freealign(ptr) ::free(ptr)
#   define SLR_alignof(T) alignof(T)
#elif defined(SLR_Platform_Linux)
#   define SLR_memalign(size, alignment) SLRAssert_NotImplemented
#   define SLR_freealign(ptr) SLRAssert_NotImplemented
#endif

// For getcwd
#if defined(SLR_Platform_Windows_MSVC)
#   define SLR_getcwd(size, buf) GetCurrentDirectory(size, buf)
#elif defined(SLR_Platform_OS_X) || defined(SLR_Platform_OpenBSD) || defined(SLR_Platform_Linux)
#   include <unistd.h>
#   define SLR_getcwd(size, buf) getcwd(buf, size)
#endif



// ----------------------------------------------------------------
// JP: よく使用される基礎的な関数の定義。
// EN: define fundamental functions often used.

template <typename T, size_t size>
constexpr size_t lengthof(const T (&array)[size]) {
    return size;
}

namespace std {
    template <typename T>
    inline T clamp(const T &v, const T &min, const T &max) {
        return std::min(max, std::max(min, v));
    }
    
    template <typename T, typename ...Args>
    inline std::array<T, sizeof...(Args)> make_array(Args &&...args) {
        return std::array<T, sizeof...(Args)>{ std::forward<Args>(args)... };
    }
}

template <typename T>
bool realEq(T a, T b, T epsilon) {
    bool forAbsolute = std::fabs(a - b) < epsilon;
    bool forRelative = std::fabs(a - b) < epsilon * std::fmax(std::fabs(a), std::fabs(b));
    return forAbsolute || forRelative;
}
template <typename T>
bool realGE(T a, T b, T epsilon) { return a > b || realEq(a - b, (T)0, epsilon); }
template <typename T>
bool realLE(T a, T b, T epsilon) { return a < b || realEq(a - b, (T)0, epsilon); }

inline uint32_t prevPowerOf2(uint32_t x) {
    x = x | (x >> 1);
    x = x | (x >> 2);
    x = x | (x >> 4);
    x = x | (x >> 8);
    x = x | (x >> 16);
    return x - (x >> 1);
}

template <typename T, typename ...ArgTypes>
std::shared_ptr<T> createShared(ArgTypes&&... args) {
    return std::shared_ptr<T>(new T(std::forward<ArgTypes>(args)...));
}

template <typename T, typename ...ArgTypes>
std::unique_ptr<T> createUnique(ArgTypes&&... args) {
    return std::unique_ptr<T>(new T(std::forward<ArgTypes>(args)...));
}

// END: define fundamental functions often used.
// ----------------------------------------------------------------



#define SLR_Color_System_CIE_1931_2deg  0
#define SLR_Color_System_CIE_1964_10deg 1
#define SLR_Color_System_CIE_2012_2deg  2
#define SLR_Color_System_CIE_2012_10deg 3



// ----------------------------------------------------------------
// JP: 機能スイッチやパラメター
// EN: Feature Switches and Parameters

#define SLR_Color_System_is_based_on SLR_Color_System_CIE_1931_2deg
#define SLR_Use_Spectral_Representation

// END: Feature Switches and Parameters
// ----------------------------------------------------------------



#if SLR_Color_System_is_based_on == SLR_Color_System_CIE_1931_2deg
#   define xbarReferenceValues xbar_CIE1931_2deg
#   define ybarReferenceValues ybar_CIE1931_2deg
#   define zbarReferenceValues zbar_CIE1931_2deg
#elif SLR_Color_System_is_based_on == SLR_Color_System_CIE_1964_10deg
#   define xbarReferenceValues xbar_CIE1964_10deg
#   define ybarReferenceValues ybar_CIE1964_10deg
#   define zbarReferenceValues zbar_CIE1964_10deg
#elif SLR_Color_System_is_based_on == SLR_Color_System_CIE_2012_2deg
#   define xbarReferenceValues xbar_CIE2012_2deg
#   define ybarReferenceValues ybar_CIE2012_2deg
#   define zbarReferenceValues zbar_CIE2012_2deg
#elif SLR_Color_System_is_based_on == SLR_Color_System_CIE_2012_10deg
#   define xbarReferenceValues xbar_CIE2012_10deg
#   define ybarReferenceValues ybar_CIE2012_10deg
#   define zbarReferenceValues zbar_CIE2012_10deg
#endif

#endif /* __SLR_defines__ */
