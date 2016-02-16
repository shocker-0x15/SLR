//
//  Vector3.h
//
//  Created by 渡部 心 on 11/07/16.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __SLR_Vector3_h__
#define __SLR_Vector3_h__

#include "../defines.h"
#include "../references.h"
#include <limits>

namespace SLR {
    template <typename RealType>
    struct SLR_API Vector3Template {
        RealType x, y, z;
        
        Vector3Template(RealType v = 0.0f) : x(v), y(v), z(v) { }
        CONSTEXPR_CONSTRUCTOR Vector3Template(RealType xx, RealType yy, RealType zz) : x(xx), y(yy), z(zz) { }
        
        Vector3Template operator+() const { return *this; }
        Vector3Template operator-() const { return Vector3Template(-x, -y, -z); }
        
        Vector3Template operator+(const Vector3Template &v) const { return Vector3Template(x + v.x, y + v.y, z + v.z); }
        Vector3Template operator-(const Vector3Template &v) const { return Vector3Template(x - v.x, y - v.y, z - v.z); }
        Vector3Template operator*(const Vector3Template &v) const { return Vector3Template(x * v.x, y * v.y, z * v.z); }
        Vector3Template operator/(const Vector3Template &v) const { return Vector3Template(x / v.x, y / v.y, z / v.z); }
        Vector3Template operator*(RealType s) const { return Vector3Template(x * s, y * s, z * s); }
        Vector3Template operator/(RealType s) const { RealType r = 1.0f / s; return Vector3Template(x * r, y * r, z * r); }
        friend inline Vector3Template operator*(RealType s, const Vector3Template &v) { return Vector3Template(s * v.x, s * v.y, s * v.z); }
        
        Vector3Template &operator+=(const Vector3Template &v) { x += v.x; y += v.y; z += v.z; return *this; }
        Vector3Template &operator-=(const Vector3Template &v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
        Vector3Template &operator*=(RealType s) { x *= s; y *= s; z *= s; return *this; }
        Vector3Template &operator/=(RealType s) { RealType r = 1.0f / s; x *= r; y *= r; z *= r; return *this; }
        
        bool operator==(const Vector3Template &v) const { return x == v.x && y == v.y && z == v.z; }
        bool operator!=(const Vector3Template &v) const { return x != v.x || y != v.y || z != v.z; }
        
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
        Vector3Template& normalize() {
            RealType length = std::sqrt(x * x + y * y + z * z);
            return *this /= length;
        }
        Vector3Template reciprocal() const { return Vector3Template(1.0f / x, 1.0f / y, 1.0f / z); }
        void makeCoordinateSystem(Vector3Template<RealType>* vx, Vector3Template<RealType>* vy) const {
            if (std::fabs(x) > std::fabs(y)) {
                float invLen = 1.0f / std::sqrt(x * x + z * z);
                *vx = Vector3Template<RealType>(-z * invLen, 0.0f, x * invLen);
            }
            else {
                float invLen = 1.0f / std::sqrt(y * y + z * z);
                *vx = Vector3Template<RealType>(0.0f, z * invLen, -y * invLen);
            }
            *vy = cross(*this, *vx);
        }
        static Vector3Template fromPolarYUp(RealType phi, RealType theta) {
            return Vector3Template(-std::sin(phi) * std::sin(theta), std::cos(theta), std::cos(phi) * std::sin(theta));
        }
        void toPolarYUp(RealType* theta, RealType* phi) const {
            *theta = std::acos(std::clamp(y, (RealType)-1.0, (RealType)1.0));
            *phi = std::fmod((RealType)(std::atan2(-x, z) + 2 * M_PI), (RealType)(2 * M_PI));
        }
        
        RealType maxValue() const { return fmaxf(x, fmaxf(y, z)); }
        RealType minValue() const { return fminf(x, fminf(y, z)); }
        bool hasNaN() const { using std::isnan; return isnan(x) || isnan(y) || isnan(z); }
        bool hasInf() const { using std::isinf; return isinf(x) || isinf(y) || isinf(z); }
        
        void print() const { printf("(%f, %f, %f)\n", x, y, z); }
        
        static const Vector3Template Zero;
        static const Vector3Template One;
        static const Vector3Template Inf;
        static const Vector3Template NaN;
        static const Vector3Template Ex;
        static const Vector3Template Ey;
        static const Vector3Template Ez;
    };
    template <typename RealType>
    const Vector3Template<RealType> Vector3Template<RealType>::Zero = Vector3Template<RealType>(0);
    template <typename RealType>
    const Vector3Template<RealType> Vector3Template<RealType>::One = Vector3Template<RealType>(1);
    template <typename RealType>
    const Vector3Template<RealType> Vector3Template<RealType>::Inf = Vector3Template<RealType>(std::numeric_limits<RealType>::infinity());
    template <typename RealType>
    const Vector3Template<RealType> Vector3Template<RealType>::NaN = Vector3Template<RealType>(std::numeric_limits<RealType>::quiet_NaN());
    template <typename RealType>
    const Vector3Template<RealType> Vector3Template<RealType>::Ex = Vector3Template<RealType>(1, 0, 0);
    template <typename RealType>
    const Vector3Template<RealType> Vector3Template<RealType>::Ey = Vector3Template<RealType>(0, 1, 0);
    template <typename RealType>
    const Vector3Template<RealType> Vector3Template<RealType>::Ez = Vector3Template<RealType>(0, 0, 1);
    
    
    template <typename RealType>
    inline Vector3Template<RealType> normalize(const Vector3Template<RealType> &v) {
        RealType l = v.length();
        return v / l;
    }
    
    template <typename RealType>
    inline RealType dot(const Vector3Template<RealType> &vec1, const Vector3Template<RealType> &vec2) {
        return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
    }
    
    template <typename RealType>
    inline Vector3Template<RealType> cross(const Vector3Template<RealType> &vec1, const Vector3Template<RealType> &vec2) {
        return Vector3Template<RealType>(vec1.y * vec2.z - vec1.z * vec2.y,
                                         vec1.z * vec2.x - vec1.x * vec2.z,
                                         vec1.x * vec2.y - vec1.y * vec2.x);
    }
    
    template <typename RealType>
    inline RealType absDot(const Vector3Template<RealType> &vec1, const Vector3Template<RealType> &vec2) {
        return std::abs(vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z);
    }
    
    template <typename RealType>
    inline Vector3Template<RealType> min(const Vector3Template<RealType> &vec1, const Vector3Template<RealType> &vec2) {
        using std::fmin;
        return Vector3Template<RealType>(fmin(vec1.x, vec2.x), fmin(vec1.y, vec2.y), fmin(vec1.z, vec2.z));
    }
    
    template <typename RealType>
    inline Vector3Template<RealType> max(const Vector3Template<RealType> &vec1, const Vector3Template<RealType> &vec2) {
        using std::fmax;
        return Vector3Template<RealType>(fmax(vec1.x, vec2.x), fmax(vec1.y, vec2.y), fmax(vec1.z, vec2.z));
    }
    
    template <typename RealType>
    inline Vector3Template<RealType> halfVector(const Vector3Template<RealType> &vec1, const Vector3Template<RealType> &vec2) {
        return normalize(vec1 + vec2);
    }
}

#endif
