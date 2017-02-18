//
//  Normal3D.h
//
//  Created by 渡部 心 on 2015/05/04.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_Normal3__
#define __SLR_Normal3__

#include "../defines.h"
#include "../declarations.h"
#include "Vector3D.h"

namespace SLR {
    template <typename RealType>
    struct SLR_API Normal3DTemplate {
        RealType x, y, z;
        
        Normal3DTemplate(RealType v = 0.0f) : x(v), y(v), z(v) { }
        CONSTEXPR_CONSTRUCTOR Normal3DTemplate(RealType xx, RealType yy, RealType zz) : x(xx), y(yy), z(zz) { }
        Normal3DTemplate(const Vector3DTemplate<RealType> &v) : x(v.x), y(v.y), z(v.z) { }
        
        operator Vector3DTemplate<RealType>() const {
            return Vector3DTemplate<RealType>(x, y, z);
        }
        
        Normal3DTemplate operator+() const { return *this; }
        Normal3DTemplate operator-() const { return Normal3DTemplate(-x, -y, -z); }
        
        Vector3DTemplate<RealType> operator+(const Vector3DTemplate<RealType> &v) const { return Vector3DTemplate<RealType>(x + v.x, y + v.y, z + v.z); }
        Vector3DTemplate<RealType> operator-(const Vector3DTemplate<RealType> &v) const { return Vector3DTemplate<RealType>(x - v.x, y - v.y, z - v.z); }
        Vector3DTemplate<RealType> operator+(const Normal3DTemplate &n) const { return Vector3DTemplate<RealType>(x + n.x, y + n.y, z + n.z); }
        Vector3DTemplate<RealType> operator-(const Normal3DTemplate &n) const { return Vector3DTemplate<RealType>(x - n.x, y - n.y, z - n.z); }
        Vector3DTemplate<RealType> operator*(RealType s) const { return Vector3DTemplate<RealType>(x * s, y * s, z * s); }
        Vector3DTemplate<RealType> operator/(RealType s) const { RealType r = 1.0f / s; return Vector3DTemplate<RealType>(x * r, y * r, z * r); }
        friend inline Vector3DTemplate<RealType> operator*(RealType s, const Normal3DTemplate &n) { return Vector3DTemplate<RealType>(s * n.x, s * n.y, s * n.z); }
        
        bool operator==(const Normal3DTemplate &p) const { return x == p.x && y == p.y && z == p.z; }
        bool operator!=(const Normal3DTemplate &p) const { return x != p.x || y != p.y || z != p.z; }
        
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
        Normal3DTemplate& normalize() {
            RealType rcpLength = 1.0 / std::sqrt(x * x + y * y + z * z);
            x *= rcpLength;
            y *= rcpLength;
            z *= rcpLength;
            return *this;
        }
        void makeCoordinateSystem(Vector3DTemplate<RealType>* tangent, Vector3DTemplate<RealType>* bitangent) const {
            if (std::fabs(x) > std::fabs(y)) {
                float invLen = 1.0f / std::sqrt(x * x + z * z);
                *tangent = Vector3DTemplate<RealType>(-z * invLen, 0.0f, x * invLen);
            }
            else {
                float invLen = 1.0f / std::sqrt(y * y + z * z);
                *tangent = Vector3DTemplate<RealType>(0.0f, z * invLen, -y * invLen);
            }
            *bitangent = cross(*this, *tangent);
        }
        
        RealType maxValue() const { using std::fmax; return fmax(x, fmax(y, z)); }
        RealType minValue() const { using std::fmin; return fmin(x, fmin(y, z)); }
        bool hasNaN() const { using std::isnan; return isnan(x) || isnan(y) || isnan(z); }
        bool hasInf() const { using std::isinf; return isinf(x) || isinf(y) || isinf(z); }
        
        std::string toString() const { char str[256]; sprintf(str, "(%g, %g, %g)", x, y, z); return str; }
        
        static const Normal3DTemplate Zero;
        static const Normal3DTemplate One;
        static const Normal3DTemplate Inf;
        static const Normal3DTemplate NaN;
    };
    template <typename RealType>
    const Normal3DTemplate<RealType> Normal3DTemplate<RealType>::Zero = Normal3DTemplate<RealType>(0);
    template <typename RealType>
    const Normal3DTemplate<RealType> Normal3DTemplate<RealType>::One = Normal3DTemplate<RealType>(1);
    template <typename RealType>
    const Normal3DTemplate<RealType> Normal3DTemplate<RealType>::Inf = Normal3DTemplate<RealType>(std::numeric_limits<RealType>::infinity());
    template <typename RealType>
    const Normal3DTemplate<RealType> Normal3DTemplate<RealType>::NaN = Normal3DTemplate<RealType>(std::numeric_limits<RealType>::quiet_NaN());
    
    
    template <typename RealType>
    inline Normal3DTemplate<RealType> normalize(const Normal3DTemplate<RealType> &n) {
        RealType l = n.length();
        return n / l;
    }
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> cross(const Vector3DTemplate<RealType> &vec, const Normal3DTemplate<RealType> &norm) {
        return Vector3DTemplate<RealType>(vec.y * norm.z - vec.z * norm.y,
                                         vec.z * norm.x - vec.x * norm.z,
                                         vec.x * norm.y - vec.y * norm.x);
    }
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> cross(const Normal3DTemplate<RealType> &norm, const Vector3DTemplate<RealType> &vec) {
        return Vector3DTemplate<RealType>(norm.y * vec.z - norm.z * vec.y,
                                         norm.z * vec.x - norm.x * vec.z,
                                         norm.x * vec.y - norm.y * vec.x);
    }
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> cross(const Normal3DTemplate<RealType> &n1, const Normal3DTemplate<RealType> &n2) {
        return Vector3DTemplate<RealType>(n1.y * n2.z - n1.z * n2.y,
                                         n1.z * n2.x - n1.x * n2.z,
                                         n1.x * n2.y - n1.y * n2.x);
    }
    
    template <typename RealType>
    inline RealType dot(const Vector3DTemplate<RealType> &vec, const Normal3DTemplate<RealType> &norm) {
        return vec.x * norm.x + vec.y * norm.y + vec.z * norm.z;
    }
    
    template <typename RealType>
    inline RealType dot(const Normal3DTemplate<RealType> &norm, const Vector3DTemplate<RealType> &vec) {
        return vec.x * norm.x + vec.y * norm.y + vec.z * norm.z;
    }
    
    template <typename RealType>
    inline RealType absDot(const Normal3DTemplate<RealType> &norm1, const Normal3DTemplate<RealType> &norm2) {
        return std::abs(norm1.x * norm2.x + norm1.y * norm2.y + norm1.z * norm2.z);
    }
    
    template <typename RealType>
    inline RealType absDot(const Vector3DTemplate<RealType> &vec, const Normal3DTemplate<RealType> &norm) {
        return std::abs(vec.x * norm.x + vec.y * norm.y + vec.z * norm.z);
    }
    
    template <typename RealType>
    inline RealType absDot(const Normal3DTemplate<RealType> &norm, const Vector3DTemplate<RealType> &vec) {
        return std::abs(vec.x * norm.x + vec.y * norm.y + vec.z * norm.z);
    }
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> min(const Normal3DTemplate<RealType> &p1, const Normal3DTemplate<RealType> &p2) {
        using std::fmin;
        return Vector3DTemplate<RealType>(fmin(p1.x, p2.x), fmin(p1.y, p2.y), fmin(p1.z, p2.z));
    }
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> max(const Normal3DTemplate<RealType> &p1, const Normal3DTemplate<RealType> &p2) {
        using std::fmax;
        return Vector3DTemplate<RealType>(fmax(p1.x, p2.x), fmax(p1.y, p2.y), fmax(p1.z, p2.z));
    }
}

#endif /* __SLR_Normal3__ */
