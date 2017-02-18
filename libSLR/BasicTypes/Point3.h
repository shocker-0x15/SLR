//
//  Point3.h
//
//  Created by 渡部 心 on 12/08/29.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __SLR_Point3__
#define __SLR_Point3__

#include "../defines.h"
#include "../references.h"
#include "Vector3.h"

namespace SLR {
    template <typename RealType>
    struct SLR_API Point3Template {
        RealType x, y, z;
        
        Point3Template(RealType v = 0.0f) : x(v), y(v), z(v) { }
        CONSTEXPR_CONSTRUCTOR Point3Template(RealType xx, RealType yy, RealType zz) : x(xx), y(yy), z(zz) { }
        Point3Template(const Vector3Template<RealType> &v) : x(v.x), y(v.y), z(v.z) { }
        
        operator Vector3Template<RealType>() const {
            return Vector3Template<RealType>(x, y, z);
        }
        
        Point3Template operator+() const { return *this; }
        Point3Template operator-() const { return Point3Template(-x, -y, -z); }
        
        Point3Template operator+(const Vector3Template<RealType> &v) const { return Point3Template(x + v.x, y + v.y, z + v.z); }
        Point3Template operator-(const Vector3Template<RealType> &v) const { return Point3Template(x - v.x, y - v.y, z - v.z); }
        Point3Template<RealType> operator+(const Point3Template &p) const { return Point3Template<RealType>(x + p.x, y + p.y, z + p.z); }
        Vector3Template<RealType> operator-(const Point3Template &p) const { return Vector3Template<RealType>(x - p.x, y - p.y, z - p.z); }
        Point3Template operator*(RealType s) const { return Point3Template(x * s, y * s, z * s); }
        Point3Template operator/(RealType s) const { RealType r = 1.0f / s; return Point3Template(x * r, y * r, z * r); }
        friend inline Point3Template operator+(const Vector3Template<RealType> &v, const Point3Template &p) { return Point3Template(p.x + v.x, p.y + v.y, p.z + v.z); }
        friend inline Point3Template operator*(RealType s, const Point3Template &p) { return Point3Template(s * p.x, s * p.y, s * p.z); }
        
        Point3Template &operator+=(const Vector3Template<RealType> &v) { x += v.x; y += v.y; z += v.z; return *this; }
        Point3Template &operator-=(const Vector3Template<RealType> &v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
        Point3Template &operator*=(RealType s) { x *= s; y *= s; z *= s; return *this; }
        Point3Template &operator/=(RealType s) { RealType r = 1.0f / s; x *= r; y *= r; z *= r; return *this; }
        
        bool operator==(const Point3Template &p) const { return x == p.x && y == p.y && z == p.z; }
        bool operator!=(const Point3Template &p) const { return x != p.x || y != p.y || z != p.z; }
        
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
        
        static const Point3Template Zero;
        static const Point3Template One;
        static const Point3Template Inf;
        static const Point3Template NaN;
    };
    template <typename RealType>
    const Point3Template<RealType> Point3Template<RealType>::Zero = Point3Template<RealType>(0);
    template <typename RealType>
    const Point3Template<RealType> Point3Template<RealType>::One = Point3Template<RealType>(1);
    template <typename RealType>
    const Point3Template<RealType> Point3Template<RealType>::Inf = Point3Template<RealType>(std::numeric_limits<RealType>::infinity());
    template <typename RealType>
    const Point3Template<RealType> Point3Template<RealType>::NaN = Point3Template<RealType>(std::numeric_limits<RealType>::quiet_NaN());
    
    
    template <typename RealType>
    inline RealType absDot(const Point3Template<RealType> &p1, const Point3Template<RealType> &p2) {
        return std::abs(p1.x * p2.x + p1.y * p2.y + p1.z * p2.z);
    }
    
    template <typename RealType>
    inline Point3Template<RealType> min(const Point3Template<RealType> &p1, const Point3Template<RealType> &p2) {
        using std::fmin;
        return Point3Template<RealType>(fmin(p1.x, p2.x), fmin(p1.y, p2.y), fmin(p1.z, p2.z));
    }
    
    template <typename RealType>
    inline Point3Template<RealType> max(const Point3Template<RealType> &p1, const Point3Template<RealType> &p2) {
        using std::fmax;
        return Point3Template<RealType>(fmax(p1.x, p2.x), fmax(p1.y, p2.y), fmax(p1.z, p2.z));
    }
    
    template <typename RealType>
    inline RealType sqDistance(const Point3Template<RealType> &p1, const Point3Template<RealType> &p2) {
        return (p2 - p1).sqLength();
    }
    
    template <typename RealType>
    inline RealType distance(const Point3Template<RealType> &p1, const Point3Template<RealType> &p2) {
        return (p2 - p1).length();
    }
}

#endif /* __SLR_Point3__ */
