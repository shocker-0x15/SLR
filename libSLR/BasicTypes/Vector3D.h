//
//  Vector3D.h
//
//  Created by 渡部 心 on 11/07/16.
//  Copyright (c) 2014年 渡部 心. All rights reserved.
//

#ifndef __SLR_Vector3D__
#define __SLR_Vector3D__

#include "../defines.h"
#include "../declarations.h"
#include <limits>

namespace SLR {
    template <typename RealType>
    struct SLR_API Vector3DTemplate {
        RealType x, y, z;
        
        Vector3DTemplate(RealType v = 0.0f) : x(v), y(v), z(v) { }
        CONSTEXPR_CONSTRUCTOR Vector3DTemplate(RealType xx, RealType yy, RealType zz) : x(xx), y(yy), z(zz) { }
        
        Vector3DTemplate operator+() const { return *this; }
        Vector3DTemplate operator-() const { return Vector3DTemplate(-x, -y, -z); }
        
        Vector3DTemplate operator+(const Vector3DTemplate &v) const { return Vector3DTemplate(x + v.x, y + v.y, z + v.z); }
        Vector3DTemplate operator-(const Vector3DTemplate &v) const { return Vector3DTemplate(x - v.x, y - v.y, z - v.z); }
        Vector3DTemplate operator*(const Vector3DTemplate &v) const { return Vector3DTemplate(x * v.x, y * v.y, z * v.z); }
        Vector3DTemplate operator/(const Vector3DTemplate &v) const { return Vector3DTemplate(x / v.x, y / v.y, z / v.z); }
        Vector3DTemplate operator*(RealType s) const { return Vector3DTemplate(x * s, y * s, z * s); }
        Vector3DTemplate operator/(RealType s) const { RealType r = 1.0f / s; return Vector3DTemplate(x * r, y * r, z * r); }
        friend inline Vector3DTemplate operator*(RealType s, const Vector3DTemplate &v) { return Vector3DTemplate(s * v.x, s * v.y, s * v.z); }
        
        Vector3DTemplate &operator+=(const Vector3DTemplate &v) { x += v.x; y += v.y; z += v.z; return *this; }
        Vector3DTemplate &operator-=(const Vector3DTemplate &v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
        Vector3DTemplate &operator*=(RealType s) { x *= s; y *= s; z *= s; return *this; }
        Vector3DTemplate &operator/=(RealType s) { RealType r = 1.0f / s; x *= r; y *= r; z *= r; return *this; }
        
        bool operator==(const Vector3DTemplate &v) const { return x == v.x && y == v.y && z == v.z; }
        bool operator!=(const Vector3DTemplate &v) const { return x != v.x || y != v.y || z != v.z; }
        
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
        Vector3DTemplate& normalize() {
            RealType length = std::sqrt(x * x + y * y + z * z);
            return *this /= length;
        }
        Vector3DTemplate reciprocal() const { return Vector3DTemplate(1.0f / x, 1.0f / y, 1.0f / z); }
        void makeCoordinateSystem(Vector3DTemplate<RealType>* vx, Vector3DTemplate<RealType>* vy) const {
            if (std::fabs(x) > std::fabs(y)) {
                float invLen = 1.0f / std::sqrt(x * x + z * z);
                *vx = Vector3DTemplate<RealType>(-z * invLen, 0.0f, x * invLen);
            }
            else {
                float invLen = 1.0f / std::sqrt(y * y + z * z);
                *vx = Vector3DTemplate<RealType>(0.0f, z * invLen, -y * invLen);
            }
            *vy = cross(*this, *vx);
        }
        static Vector3DTemplate fromPolarZUp(RealType phi, RealType theta) {
            return Vector3DTemplate(std::sin(phi) * std::sin(theta), std::cos(phi) * std::sin(theta), std::cos(theta));
        }
        static Vector3DTemplate fromPolarYUp(RealType phi, RealType theta) {
            return Vector3DTemplate(-std::sin(phi) * std::sin(theta), std::cos(theta), std::cos(phi) * std::sin(theta));
        }
        void toPolarYUp(RealType* theta, RealType* phi) const {
            *theta = std::acos(std::clamp(y, (RealType)-1.0, (RealType)1.0));
            *phi = std::fmod((RealType)(std::atan2(-x, z) + 2 * M_PI), (RealType)(2 * M_PI));
        }
        
        RealType maxValue() const { return fmaxf(x, fmaxf(y, z)); }
        RealType minValue() const { return fminf(x, fminf(y, z)); }
        bool hasNaN() const { using std::isnan; return isnan(x) || isnan(y) || isnan(z); }
        bool hasInf() const { using std::isinf; return isinf(x) || isinf(y) || isinf(z); }
        
        std::string toString() const { char str[256]; sprintf(str, "(%g, %g, %g)", x, y, z); return str; }
        
        static const Vector3DTemplate Zero;
        static const Vector3DTemplate One;
        static const Vector3DTemplate Inf;
        static const Vector3DTemplate NaN;
        static const Vector3DTemplate Ex;
        static const Vector3DTemplate Ey;
        static const Vector3DTemplate Ez;
    };
    template <typename RealType>
    const Vector3DTemplate<RealType> Vector3DTemplate<RealType>::Zero = Vector3DTemplate<RealType>(0);
    template <typename RealType>
    const Vector3DTemplate<RealType> Vector3DTemplate<RealType>::One = Vector3DTemplate<RealType>(1);
    template <typename RealType>
    const Vector3DTemplate<RealType> Vector3DTemplate<RealType>::Inf = Vector3DTemplate<RealType>(std::numeric_limits<RealType>::infinity());
    template <typename RealType>
    const Vector3DTemplate<RealType> Vector3DTemplate<RealType>::NaN = Vector3DTemplate<RealType>(std::numeric_limits<RealType>::quiet_NaN());
    template <typename RealType>
    const Vector3DTemplate<RealType> Vector3DTemplate<RealType>::Ex = Vector3DTemplate<RealType>(1, 0, 0);
    template <typename RealType>
    const Vector3DTemplate<RealType> Vector3DTemplate<RealType>::Ey = Vector3DTemplate<RealType>(0, 1, 0);
    template <typename RealType>
    const Vector3DTemplate<RealType> Vector3DTemplate<RealType>::Ez = Vector3DTemplate<RealType>(0, 0, 1);
    
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> normalize(const Vector3DTemplate<RealType> &v) {
        RealType l = v.length();
        return v / l;
    }
    
    template <typename RealType>
    inline RealType dot(const Vector3DTemplate<RealType> &vec1, const Vector3DTemplate<RealType> &vec2) {
        return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
    }
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> cross(const Vector3DTemplate<RealType> &vec1, const Vector3DTemplate<RealType> &vec2) {
        return Vector3DTemplate<RealType>(vec1.y * vec2.z - vec1.z * vec2.y,
                                         vec1.z * vec2.x - vec1.x * vec2.z,
                                         vec1.x * vec2.y - vec1.y * vec2.x);
    }
    
    template <typename RealType>
    inline RealType absDot(const Vector3DTemplate<RealType> &vec1, const Vector3DTemplate<RealType> &vec2) {
        return std::abs(vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z);
    }
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> min(const Vector3DTemplate<RealType> &vec1, const Vector3DTemplate<RealType> &vec2) {
        using std::fmin;
        return Vector3DTemplate<RealType>(fmin(vec1.x, vec2.x), fmin(vec1.y, vec2.y), fmin(vec1.z, vec2.z));
    }
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> max(const Vector3DTemplate<RealType> &vec1, const Vector3DTemplate<RealType> &vec2) {
        using std::fmax;
        return Vector3DTemplate<RealType>(fmax(vec1.x, vec2.x), fmax(vec1.y, vec2.y), fmax(vec1.z, vec2.z));
    }
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> halfVector(const Vector3DTemplate<RealType> &vec1, const Vector3DTemplate<RealType> &vec2) {
        return normalize(vec1 + vec2);
    }
}

#endif /* __SLR_Vector3D__ */
