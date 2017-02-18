//
//  Vector4D.h
//
//  Created by 渡部 心 on 11/08/22.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __SLR_Vector4D__
#define __SLR_Vector4D__

#include "../defines.h"
#include "../declarations.h"
#include "Vector3D.h"

namespace SLR {
    template <typename RealType>
    struct SLR_API Vector4DTemplate {
    public:
        RealType x, y, z, w;
        
        Vector4DTemplate(RealType v = 0.0f) : x(v), y(v), z(v), w(v) { }
        CONSTEXPR_CONSTRUCTOR Vector4DTemplate(RealType xx, RealType yy, RealType zz, RealType ww) : x(xx), y(yy), z(zz), w(ww) { }
        Vector4DTemplate(const Vector3DTemplate<RealType> &vec3, RealType ww) : x(vec3.x), y(vec3.y), z(vec3.z), w(ww) { }
        
        Vector4DTemplate operator+() const { return *this; }
        Vector4DTemplate operator-() const { return Vector4DTemplate(-x, -y, -z, -w); }
        
        Vector4DTemplate operator+(const Vector4DTemplate &v) const { return Vector4DTemplate(x + v.x, y + v.y, z + v.z, w + v.w); }
        Vector4DTemplate operator-(const Vector4DTemplate &v) const { return Vector4DTemplate(x - v.x, y - v.y, z - v.z, w - v.w); }
        Vector4DTemplate operator*(RealType s) const { return Vector4DTemplate(x * s, y * s, z * s,  w * s); }
        Vector4DTemplate operator/(RealType s) const { RealType r = 1.0f / s; return Vector4DTemplate(x * r, y * r, z * r, w * r); }
        friend inline Vector4DTemplate operator*(RealType s, const Vector4DTemplate &v) { return Vector4DTemplate(s * v.x, s * v.y, s * v.z, s * v.w); }
        
        Vector4DTemplate &operator+=(const Vector4DTemplate &v) { x += v.x; y += v.y; z += v.z; w += v.w; return *this; }
        Vector4DTemplate &operator-=(const Vector4DTemplate &v) { x -= v.x; y -= v.y; z -= v.z; w -= v.w; return *this; }
        Vector4DTemplate &operator*=(RealType s) { x *= s; y *= s; z *= s; w *= s; return *this; }
        Vector4DTemplate &operator/=(RealType s) { RealType r = 1.0f / s; x *= r; y *= r; z *= r; w *= r; return *this; }
        
        bool operator==(const Vector4DTemplate &v) const { return x == v.x && y == v.y && z == v.z && w == v.w; }
        bool operator!=(const Vector4DTemplate &v) const { return x != v.x || y != v.y || z != v.z || w != v.w; }
        
        RealType &operator[](unsigned int index) {
            SLRAssert(index < 4, "\"index\" is out of range [0, 3].");
            return *(&x + index);
        }
        RealType operator[](unsigned int index) const {
            SLRAssert(index < 4, "\"index\" is out of range [0, 3].");
            return *(&x + index);
        }
        
        explicit operator Vector3DTemplate<RealType>() const { return Vector3DTemplate<RealType>(x, y, z); }
        
        RealType maxValue() const { using std::fmax; return fmax(fmax(x, y), fmax(z, w)); }
        RealType minValue() const { using std::fmin; return fmin(fmin(x, y), fmin(z, w)); }
        bool hasNaN() const { using std::isnan; return isnan(x) || isnan(y) || isnan(z) || isnan(w); }
        bool hasInf() const { using std::isinf; return isinf(x) || isinf(y) || isinf(z) || isinf(w); }
        
        std::string toString() const { char str[256]; sprintf(str, "(%g, %g, %g, %g)", x, y, z, w); return str; }
        
        static const Vector4DTemplate Zero;
        static const Vector4DTemplate One;
        static const Vector4DTemplate Inf;
        static const Vector4DTemplate NaN;
        static const Vector4DTemplate Ex;
        static const Vector4DTemplate Ey;
        static const Vector4DTemplate Ez;
        static const Vector4DTemplate Ew;
    };
    template <typename RealType>
    const Vector4DTemplate<RealType> Vector4DTemplate<RealType>::Zero = Vector4DTemplate<RealType>(0);
    template <typename RealType>
    const Vector4DTemplate<RealType> Vector4DTemplate<RealType>::One = Vector4DTemplate<RealType>(1);
    template <typename RealType>
    const Vector4DTemplate<RealType> Vector4DTemplate<RealType>::Inf = Vector4DTemplate<RealType>(std::numeric_limits<RealType>::infinity());
    template <typename RealType>
    const Vector4DTemplate<RealType> Vector4DTemplate<RealType>::NaN = Vector4DTemplate<RealType>(std::numeric_limits<RealType>::quiet_NaN());
    template <typename RealType>
    const Vector4DTemplate<RealType> Vector4DTemplate<RealType>::Ex = Vector4DTemplate<RealType>(1, 0, 0, 0);
    template <typename RealType>
    const Vector4DTemplate<RealType> Vector4DTemplate<RealType>::Ey = Vector4DTemplate<RealType>(0, 1, 0, 0);
    template <typename RealType>
    const Vector4DTemplate<RealType> Vector4DTemplate<RealType>::Ez = Vector4DTemplate<RealType>(0, 0, 1, 0);
    template <typename RealType>
    const Vector4DTemplate<RealType> Vector4DTemplate<RealType>::Ew = Vector4DTemplate<RealType>(0, 0, 0, 1);
    
    
    template <typename RealType>
    inline RealType dot(const Vector4DTemplate<RealType> &vec1, const Vector4DTemplate<RealType> &vec2) {
        return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z + vec1.w * vec2.w;
    }
    
    template <typename RealType>
    inline Vector4DTemplate<RealType> min(const Vector4DTemplate<RealType> &vec1, const Vector4DTemplate<RealType> &vec2) {
        using std::fmin;
        return Vector4DTemplate<RealType>(fmin(vec1.x, vec2.x), fmin(vec1.y, vec2.y), fmin(vec1.z, vec2.z), fmin(vec1.w, vec2.w));
    }
    
    template <typename RealType>
    inline Vector4DTemplate<RealType> max(const Vector4DTemplate<RealType> &vec1, const Vector4DTemplate<RealType> &vec2) {
        using std::fmax;
        return Vector4DTemplate<RealType>(fmax(vec1.x, vec2.x), fmax(vec1.y, vec2.y), fmax(vec1.z, vec2.z), fmax(vec1.w, vec2.w));
    }    
}

#endif /* __SLR_Vector4D__ */
