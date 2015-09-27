//
//  Vector4.h
//
//  Created by 渡部 心 on 11/08/22.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __SLR_Vector4_h__
#define __SLR_Vector4_h__

#include "../defines.h"
#include "../references.h"
#include "Vector3.h"

namespace SLR {
    template <typename RealType>
    struct Vector4Template {
    public:
        RealType x, y, z, w;
        
        Vector4Template(float v = 0.0f) : x(v), y(v), z(v), w(v) { };
        constexpr Vector4Template(RealType xx, RealType yy, RealType zz, RealType ww) : x(xx), y(yy), z(zz), w(ww) { };
        Vector4Template(const Vector3Template<RealType> &vec3, RealType ww) : x(vec3.x), y(vec3.y), z(vec3.z), w(ww) { };
        
        Vector4Template operator+() const { return *this; };
        Vector4Template operator-() const { return Vector4Template(-x, -y, -z, -w); };
        
        Vector4Template operator+(const Vector4Template &v) const { return Vector4Template(x + v.x, y + v.y, z + v.z, w + v.w); }
        Vector4Template operator-(const Vector4Template &v) const { return Vector4Template(x - v.x, y - v.y, z - v.z, w - v.w); };
        Vector4Template operator*(RealType s) const { return Vector4Template(x * s, y * s, z * s,  w * s); };
        Vector4Template operator/(RealType s) const { RealType r = 1.0f / s; return Vector4Template(x * r, y * r, z * r, w * r); };
        friend inline Vector4Template operator*(RealType s, const Vector4Template &v) { return Vector4Template(s * v.x, s * v.y, s * v.z, s * v.w); };
        
        Vector4Template &operator+=(const Vector4Template &v) { x += v.x; y += v.y; z += v.z; w += v.w; return *this; };
        Vector4Template &operator-=(const Vector4Template &v) { x -= v.x; y -= v.y; z -= v.z; w -= v.w; return *this; };
        Vector4Template &operator*=(RealType s) { x *= s; y *= s; z *= s; w *= s; return *this; };
        Vector4Template &operator/=(RealType s) { RealType r = 1.0f / s; x *= r; y *= r; z *= r; w *= r; return *this; };
        
        bool operator==(const Vector4Template &v) const { return x == v.x && y == v.y && z == v.z && w == v.w; };
        bool operator!=(const Vector4Template &v) const { return x != v.x || y != v.y || z != v.z || w != v.w; };
        
        RealType &operator[](unsigned int index) {
            SLRAssert(index < 4, "\"index\" is out of range [0, 3].");
            return *(&x + index);
        };
        RealType operator[](unsigned int index) const {
            SLRAssert(index < 4, "\"index\" is out of range [0, 3].");
            return *(&x + index);
        };
        
        explicit operator Vector3Template<RealType>() const { return Vector3Template<RealType>(x, y, z); };
        
        RealType maxValue() const { return fmaxf(fmaxf(x, y), fmaxf(z, w)); };
        RealType minValue() const { return fminf(fminf(x, y), fminf(z, w)); };
        bool hasNaN() const { using std::isnan; return isnan(x) || isnan(y) || isnan(z) || isnan(w); };
        bool hasInf() const { using std::isinf; return isinf(x) || isinf(y) || isinf(z) || isinf(w); };
        
        static const Vector4Template Zero;
        static const Vector4Template One;
        static const Vector4Template Inf;
        static const Vector4Template NaN;
        static const Vector4Template Ex;
        static const Vector4Template Ey;
        static const Vector4Template Ez;
        static const Vector4Template Ew;
    };
    
    template <typename RealType>
    const Vector4Template<RealType> Vector4Template<RealType>::Zero = Vector4Template<RealType>(0);
    template <typename RealType>
    const Vector4Template<RealType> Vector4Template<RealType>::One = Vector4Template<RealType>(1);
    template <typename RealType>
    const Vector4Template<RealType> Vector4Template<RealType>::Inf = Vector4Template<RealType>(std::numeric_limits<RealType>::infinity());
    template <typename RealType>
    const Vector4Template<RealType> Vector4Template<RealType>::NaN = Vector4Template<RealType>(std::numeric_limits<RealType>::quiet_NaN());
    template <typename RealType>
    const Vector4Template<RealType> Vector4Template<RealType>::Ex = Vector4Template<RealType>(1, 0, 0, 0);
    template <typename RealType>
    const Vector4Template<RealType> Vector4Template<RealType>::Ey = Vector4Template<RealType>(0, 1, 0, 0);
    template <typename RealType>
    const Vector4Template<RealType> Vector4Template<RealType>::Ez = Vector4Template<RealType>(0, 0, 1, 0);
    template <typename RealType>
    const Vector4Template<RealType> Vector4Template<RealType>::Ew = Vector4Template<RealType>(0, 0, 0, 1);
    
    
    template <typename RealType>
    inline RealType dot(const Vector4Template<RealType> &vec1, const Vector4Template<RealType> &vec2) {
        return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z + vec1.w * vec2.w;
    }
    
    template <typename RealType>
    inline Vector4Template<RealType> min(const Vector4Template<RealType> &vec1, const Vector4Template<RealType> &vec2) {
        using std::fmin;
        return Vector4Template<RealType>(fmin(vec1.x, vec2.x), fmin(vec1.y, vec2.y), fmin(vec1.z, vec2.z), fmin(vec1.w, vec2.w));
    }
    
    template <typename RealType>
    inline Vector4Template<RealType> max(const Vector4Template<RealType> &vec1, const Vector4Template<RealType> &vec2) {
        using std::fmax;
        return Vector4Template<RealType>(fmax(vec1.x, vec2.x), fmax(vec1.y, vec2.y), fmax(vec1.z, vec2.z), fmax(vec1.w, vec2.w));
    }    
}

#endif
