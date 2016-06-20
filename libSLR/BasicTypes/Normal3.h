//
//  Normal3.h
//
//  Created by 渡部 心 on 2015/05/04.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_Normal3_h__
#define __SLR_Normal3_h__

#include "../defines.h"
#include "../references.h"
#include "Vector3.h"

namespace SLR {
    template <typename RealType>
    struct SLR_API Normal3Template {
        RealType x, y, z;
        
        Normal3Template(RealType v = 0.0f) : x(v), y(v), z(v) { }
        CONSTEXPR_CONSTRUCTOR Normal3Template(RealType xx, RealType yy, RealType zz) : x(xx), y(yy), z(zz) { }
        Normal3Template(const Vector3Template<RealType> &v) : x(v.x), y(v.y), z(v.z) { }
        
        operator Vector3Template<RealType>() const {
            return Vector3Template<RealType>(x, y, z);
        }
        
        Normal3Template operator+() const { return *this; }
        Normal3Template operator-() const { return Normal3Template(-x, -y, -z); }
        
        Vector3Template<RealType> operator+(const Vector3Template<RealType> &v) const { return Vector3Template<RealType>(x + v.x, y + v.y, z + v.z); }
        Vector3Template<RealType> operator-(const Vector3Template<RealType> &v) const { return Vector3Template<RealType>(x - v.x, y - v.y, z - v.z); }
        Vector3Template<RealType> operator+(const Normal3Template &n) const { return Vector3Template<RealType>(x + n.x, y + n.y, z + n.z); }
        Vector3Template<RealType> operator-(const Normal3Template &n) const { return Vector3Template<RealType>(x - n.x, y - n.y, z - n.z); }
        Vector3Template<RealType> operator*(RealType s) const { return Vector3Template<RealType>(x * s, y * s, z * s); }
        Vector3Template<RealType> operator/(RealType s) const { RealType r = 1.0f / s; return Vector3Template<RealType>(x * r, y * r, z * r); }
        friend inline Vector3Template<RealType> operator*(RealType s, const Normal3Template &n) { return Vector3Template<RealType>(s * n.x, s * n.y, s * n.z); }
        
        bool operator==(const Normal3Template &p) const { return x == p.x && y == p.y && z == p.z; }
        bool operator!=(const Normal3Template &p) const { return x != p.x || y != p.y || z != p.z; }
        
        RealType &operator[](unsigned int index) {
            SLRAssert(index < 3, "\"index\" is out of range [0, 2].");
            return *(&x + index);
        }
        RealType operator[](unsigned int index) const {
            SLRAssert(index < 3, "\"index\" is out of range [0, 2].");
            return *(&x + index);
        }
        
        RealType length() const { return std::sqrt(x * x + y * y + z * z); }
        RealType sqLength() const { return x * x + y * y + z * z; }
        Normal3Template& normalize() {
            RealType rcpLength = 1.0 / std::sqrt(x * x + y * y + z * z);
            x *= rcpLength;
            y *= rcpLength;
            z *= rcpLength;
            return *this;
        }
        void makeCoordinateSystem(Vector3Template<RealType>* tangent, Vector3Template<RealType>* bitangent) const {
            if (std::fabs(x) > std::fabs(y)) {
                float invLen = 1.0f / std::sqrt(x * x + z * z);
                *tangent = Vector3Template<RealType>(-z * invLen, 0.0f, x * invLen);
            }
            else {
                float invLen = 1.0f / std::sqrt(y * y + z * z);
                *tangent = Vector3Template<RealType>(0.0f, z * invLen, -y * invLen);
            }
            *bitangent = cross(*this, *tangent);
        }
        
        RealType maxValue() const { using std::fmax; return fmax(x, fmax(y, z)); }
        RealType minValue() const { using std::fmin; return fmin(x, fmin(y, z)); }
        bool hasNaN() const { using std::isnan; return isnan(x) || isnan(y) || isnan(z); }
        bool hasInf() const { using std::isinf; return isinf(x) || isinf(y) || isinf(z); }
        
        std::string toString() const { char str[256]; sprintf(str, "(%g, %g, %g)", x, y, z); return str; }
        
        static const Normal3Template Zero;
        static const Normal3Template One;
        static const Normal3Template Inf;
        static const Normal3Template NaN;
    };
    template <typename RealType>
    const Normal3Template<RealType> Normal3Template<RealType>::Zero = Normal3Template<RealType>(0);
    template <typename RealType>
    const Normal3Template<RealType> Normal3Template<RealType>::One = Normal3Template<RealType>(1);
    template <typename RealType>
    const Normal3Template<RealType> Normal3Template<RealType>::Inf = Normal3Template<RealType>(std::numeric_limits<RealType>::infinity());
    template <typename RealType>
    const Normal3Template<RealType> Normal3Template<RealType>::NaN = Normal3Template<RealType>(std::numeric_limits<RealType>::quiet_NaN());
    
    
    template <typename RealType>
    inline Normal3Template<RealType> normalize(const Normal3Template<RealType> &n) {
        RealType l = n.length();
        return n / l;
    }
    
    template <typename RealType>
    inline Vector3Template<RealType> cross(const Vector3Template<RealType> &vec, const Normal3Template<RealType> &norm) {
        return Vector3Template<RealType>(vec.y * norm.z - vec.z * norm.y,
                                         vec.z * norm.x - vec.x * norm.z,
                                         vec.x * norm.y - vec.y * norm.x);
    }
    
    template <typename RealType>
    inline Vector3Template<RealType> cross(const Normal3Template<RealType> &norm, const Vector3Template<RealType> &vec) {
        return Vector3Template<RealType>(norm.y * vec.z - norm.z * vec.y,
                                         norm.z * vec.x - norm.x * vec.z,
                                         norm.x * vec.y - norm.y * vec.x);
    }
    
    template <typename RealType>
    inline Vector3Template<RealType> cross(const Normal3Template<RealType> &n1, const Normal3Template<RealType> &n2) {
        return Vector3Template<RealType>(n1.y * n2.z - n1.z * n2.y,
                                         n1.z * n2.x - n1.x * n2.z,
                                         n1.x * n2.y - n1.y * n2.x);
    }
    
    template <typename RealType>
    inline RealType dot(const Vector3Template<RealType> &vec, const Normal3Template<RealType> &norm) {
        return vec.x * norm.x + vec.y * norm.y + vec.z * norm.z;
    }
    
    template <typename RealType>
    inline RealType dot(const Normal3Template<RealType> &norm, const Vector3Template<RealType> &vec) {
        return vec.x * norm.x + vec.y * norm.y + vec.z * norm.z;
    }
    
    template <typename RealType>
    inline RealType absDot(const Normal3Template<RealType> &norm1, const Normal3Template<RealType> &norm2) {
        return std::abs(norm1.x * norm2.x + norm1.y * norm2.y + norm1.z * norm2.z);
    }
    
    template <typename RealType>
    inline RealType absDot(const Vector3Template<RealType> &vec, const Normal3Template<RealType> &norm) {
        return std::abs(vec.x * norm.x + vec.y * norm.y + vec.z * norm.z);
    }
    
    template <typename RealType>
    inline RealType absDot(const Normal3Template<RealType> &norm, const Vector3Template<RealType> &vec) {
        return std::abs(vec.x * norm.x + vec.y * norm.y + vec.z * norm.z);
    }
    
    template <typename RealType>
    inline Vector3Template<RealType> min(const Normal3Template<RealType> &p1, const Normal3Template<RealType> &p2) {
        using std::fmin;
        return Vector3Template<RealType>(fmin(p1.x, p2.x), fmin(p1.y, p2.y), fmin(p1.z, p2.z));
    }
    
    template <typename RealType>
    inline Vector3Template<RealType> max(const Normal3Template<RealType> &p1, const Normal3Template<RealType> &p2) {
        using std::fmax;
        return Vector3Template<RealType>(fmax(p1.x, p2.x), fmax(p1.y, p2.y), fmax(p1.z, p2.z));
    }
}

#endif
