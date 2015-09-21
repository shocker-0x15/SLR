//
//  defines.h
//
//  Created by 渡部 心 on 2015/07/08.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLR_defines_h
#define SLR_defines_h

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include <vector>
#include <map>
#include <set>
#include <stack>

#include <limits>
#include <algorithm>
#include <memory>
#include <functional>

#ifdef DEBUG
#define ENABLE_ASSERT
#endif

#ifdef ENABLE_ASSERT
#    define SLRAssert(expr, fmt, ...) if (!(expr)) { printf("%s @%s: %u:\n", #expr, __FILE__, __LINE__); printf(fmt"\n", ##__VA_ARGS__); abort(); } 0
#else
#    define SLRAssert(expr, fmt, ...)
#endif

#define SLRAssert_NotDefined SLRAssert(false, "Not defined!");

// Platform defines
#if defined(_WIN32) || defined(_WIN64)
#    define SLR_Defs_Windows
#    if defined(__MINGW32__) // Defined for both 32 bit/64 bit MinGW
#        define SLR_Defs_MinGW
#    elif defined(_MSC_VER)
#        define SLR_Defs_MSVC
#    endif
#elif defined(__linux__)
#    define SLR_Defs_Linux
#elif defined(__APPLE__)
#    define SLR_Defs_OS_X
#elif defined(__OpenBSD__)
#    define SLR_Defs_OpenBSD
#endif

#define SLR_Minimum_Machine_Alignment 16
#define SLR_L1_Cacheline_Size 64

// For memalign, free, alignof
#if defined(SLR_Defs_Windows)
#    include <malloc.h>
#    define SLR_memalign(size, alignment) _aligned_malloc(size, alignment)
#    define SLR_freealign(ptr) _aligned_free(ptr)
#    define SLR_alignof(T) __alignof(T)
#elif defined(SLR_Defs_OS_X) || defined(SLR_Defs_OpenBSD)
inline void* SLR_memalign(size_t size, size_t alignment) {
    void* ptr;
    if (posix_memalign(&ptr, alignment, size))
        ptr = nullptr;
    return ptr;
}
#     define SLR_freealign(ptr) ::free(ptr)
#     define SLR_alignof(T) alignof(T)
#elif defined(SLR_Defs_Linux)
#    define SLR_memalign(size, alignment) SLRAssert_NotDefined
#     define SLR_freealign(ptr) SLRAssert_NotDefined
#endif

namespace std {
    template <typename T>
    inline T clamp(const T &v, const T &min, const T &max) {
        return std::min(max, std::max(min, v));
    }
}

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

#define Use_BSDF_Actual_Weights
#define Use_Spectral_Representation

#endif
