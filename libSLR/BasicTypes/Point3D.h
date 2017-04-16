//
//  Point3D.h
//
//  Created by 渡部 心 on 12/08/29.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __SLR_Point3D__
#define __SLR_Point3D__

#include "../defines.h"
#include "../declarations.h"
#include "Vector3D.h"

namespace SLR {
    template <typename RealType>
    struct SLR_API Point3DTemplate {
        RealType x, y, z;
        
        Point3DTemplate(RealType v = 0.0f) : x(v), y(v), z(v) { }
        CONSTEXPR_CONSTRUCTOR Point3DTemplate(RealType xx, RealType yy, RealType zz) : x(xx), y(yy), z(zz) { }
        Point3DTemplate(const Vector3DTemplate<RealType> &v) : x(v.x), y(v.y), z(v.z) { }
        
        operator Vector3DTemplate<RealType>() const {
            return Vector3DTemplate<RealType>(x, y, z);
        }
        
        Point3DTemplate operator+() const { return *this; }
        Point3DTemplate operator-() const { return Point3DTemplate(-x, -y, -z); }
        
        Point3DTemplate operator+(const Vector3DTemplate<RealType> &v) const { return Point3DTemplate(x + v.x, y + v.y, z + v.z); }
        Point3DTemplate operator-(const Vector3DTemplate<RealType> &v) const { return Point3DTemplate(x - v.x, y - v.y, z - v.z); }
        Point3DTemplate<RealType> operator+(const Point3DTemplate &p) const { return Point3DTemplate<RealType>(x + p.x, y + p.y, z + p.z); }
        Vector3DTemplate<RealType> operator-(const Point3DTemplate &p) const { return Vector3DTemplate<RealType>(x - p.x, y - p.y, z - p.z); }
        Point3DTemplate operator*(RealType s) const { return Point3DTemplate(x * s, y * s, z * s); }
        Point3DTemplate operator/(RealType s) const { RealType r = 1.0f / s; return Point3DTemplate(x * r, y * r, z * r); }
        friend inline Point3DTemplate operator+(const Vector3DTemplate<RealType> &v, const Point3DTemplate &p) { return Point3DTemplate(p.x + v.x, p.y + v.y, p.z + v.z); }
        friend inline Point3DTemplate operator*(RealType s, const Point3DTemplate &p) { return Point3DTemplate(s * p.x, s * p.y, s * p.z); }
        
        Point3DTemplate &operator+=(const Vector3DTemplate<RealType> &v) { x += v.x; y += v.y; z += v.z; return *this; }
        Point3DTemplate &operator-=(const Vector3DTemplate<RealType> &v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
        Point3DTemplate &operator*=(RealType s) { x *= s; y *= s; z *= s; return *this; }
        Point3DTemplate &operator/=(RealType s) { RealType r = 1.0f / s; x *= r; y *= r; z *= r; return *this; }
        
        bool operator==(const Point3DTemplate &p) const { return x == p.x && y == p.y && z == p.z; }
        bool operator!=(const Point3DTemplate &p) const { return x != p.x || y != p.y || z != p.z; }
        
        RealType &operator[](unsigned int index) {
            SLRAssert(index < 3, "\"index\" is out of range [0, 2].");
            return *(&x + index);
        }
        RealType operator[](unsigned int index) const {
            SLRAssert(index < 3, "\"index\" is out of range [0, 2].");
            return *(&x + index);
        }
        
        RealType maxValue() const { return std::fmax(x, std::fmax(y, z)); }
        RealType minValue() const { return std::fmin(x, std::fmin(y, z)); }
        bool hasNaN() const { using std::isnan; return isnan(x) || isnan(y) || isnan(z); }
        bool hasInf() const { using std::isinf; return isinf(x) || isinf(y) || isinf(z); }
        
        std::string toString() const { char str[256]; sprintf(str, "(%g, %g, %g)", x, y, z); return str; }
        
        static const Point3DTemplate Zero;
        static const Point3DTemplate One;
        static const Point3DTemplate Inf;
        static const Point3DTemplate NaN;
    };
    template <typename RealType>
    const Point3DTemplate<RealType> Point3DTemplate<RealType>::Zero = Point3DTemplate<RealType>(0);
    template <typename RealType>
    const Point3DTemplate<RealType> Point3DTemplate<RealType>::One = Point3DTemplate<RealType>(1);
    template <typename RealType>
    const Point3DTemplate<RealType> Point3DTemplate<RealType>::Inf = Point3DTemplate<RealType>(std::numeric_limits<RealType>::infinity());
    template <typename RealType>
    const Point3DTemplate<RealType> Point3DTemplate<RealType>::NaN = Point3DTemplate<RealType>(std::numeric_limits<RealType>::quiet_NaN());
    
    
    template <typename RealType>
    inline RealType absDot(const Point3DTemplate<RealType> &p1, const Point3DTemplate<RealType> &p2) {
        return std::abs(p1.x * p2.x + p1.y * p2.y + p1.z * p2.z);
    }
    
    template <typename RealType>
    inline Point3DTemplate<RealType> min(const Point3DTemplate<RealType> &p1, const Point3DTemplate<RealType> &p2) {
        using std::fmin;
        return Point3DTemplate<RealType>(fmin(p1.x, p2.x), fmin(p1.y, p2.y), fmin(p1.z, p2.z));
    }
    
    template <typename RealType>
    inline Point3DTemplate<RealType> max(const Point3DTemplate<RealType> &p1, const Point3DTemplate<RealType> &p2) {
        using std::fmax;
        return Point3DTemplate<RealType>(fmax(p1.x, p2.x), fmax(p1.y, p2.y), fmax(p1.z, p2.z));
    }
    
    template <typename RealType>
    inline Point3DTemplate<RealType> clamp(const Point3DTemplate<RealType> &p, const Point3DTemplate<RealType> &minP, const Point3DTemplate<RealType> &maxP) {
        return max(min(p, maxP), minP);
    }
    
    template <typename RealType>
    inline RealType sqDistance(const Point3DTemplate<RealType> &p1, const Point3DTemplate<RealType> &p2) {
        return (p2 - p1).sqLength();
    }
    
    template <typename RealType>
    inline RealType distance(const Point3DTemplate<RealType> &p1, const Point3DTemplate<RealType> &p2) {
        return (p2 - p1).length();
    }
}

#endif /* __SLR_Point3D__ */
