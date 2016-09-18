//
//  Matrix3x3.h
//
//  Created by 渡部 心 on 2016/09/08.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_Matrix3x3_h__
#define __SLR_Matrix3x3_h__

#include "../defines.h"
#include "Vector3.h"
#include "Point3.h"

namespace SLR {
    template <typename RealType>
    struct SLR_API Matrix3x3Template {
        union {
            struct { RealType m00, m10, m20; };
            Vector3Template<RealType> c0;
        };
        union {
            struct { RealType m01, m11, m21; };
            Vector3Template<RealType> c1;
        };
        union {
            struct { RealType m02, m12, m22; };
            Vector3Template<RealType> c2;
        };
        
        Matrix3x3Template() : c0(Vector3Template<RealType>::Zero), c1(Vector3Template<RealType>::Zero), c2(Vector3Template<RealType>::Zero) { }
        Matrix3x3Template(RealType array[16]) :
        m00(array[0]), m10(array[1]), m20(array[2]),
        m01(array[3]), m11(array[4]), m21(array[5]),
        m02(array[6]), m12(array[7]), m22(array[8]) { }
        CONSTEXPR_CONSTRUCTOR Matrix3x3Template(const Vector3Template<RealType> &col0, const Vector3Template<RealType> &col1, const Vector3Template<RealType> &col2) :
        c0(col0), c1(col1), c2(col2)
#ifdef SLR_Platform_Windows_MSVC
        ,// without these extra initialization, an error occurs at least with VS2015.
        m00(col0.x), m10(col0.y), m20(col0.z),
        m01(col1.x), m11(col1.y), m21(col1.z),
        m02(col2.x), m12(col2.y), m22(col2.z)
#endif
        { }
        
        Matrix3x3Template operator+() const { return *this; }
        Matrix3x3Template operator-() const { return Matrix3x3Template(-c0, -c1, -c2); }
        
        Matrix3x3Template operator+(const Matrix3x3Template &mat) const { return Matrix3x3Template(c0 + mat.c0, c1 + mat.c1, c2 + mat.c2); }
        Matrix3x3Template operator-(const Matrix3x3Template &mat) const { return Matrix3x3Template(c0 - mat.c0, c1 - mat.c1, c2 - mat.c2); }
        Matrix3x3Template operator*(const Matrix3x3Template &mat) const {
            const Vector3Template<RealType> r[] = {row(0), row(1), row(2)};
            return Matrix3x3Template(Vector3Template<RealType>(dot(r[0], mat.c0), dot(r[1], mat.c0), dot(r[2], mat.c0)),
                                     Vector3Template<RealType>(dot(r[0], mat.c1), dot(r[1], mat.c1), dot(r[2], mat.c1)),
                                     Vector3Template<RealType>(dot(r[0], mat.c2), dot(r[1], mat.c2), dot(r[2], mat.c2)));
        }
        Vector3Template<RealType> operator*(const Vector3Template<RealType> &v) const {
            return Vector3Template<RealType>(dot(row(0), v), dot(row(1), v), dot(row(2), v));
        }
        Point3Template<RealType> operator*(const Point3Template<RealType> &p) const {
            Vector3Template<RealType> ph{p.x, p.y, p.z};
            Vector3Template<RealType> pht = Vector3Template<RealType>(dot(row(0), ph), dot(row(1), ph), dot(row(2), ph));
            return Point3Template<RealType>(pht.x, pht.y, pht.z);
        }
        Matrix3x3Template operator*(RealType s) const { return Matrix3x3Template(c0 * s, c1 * s, c2 * s); }
        Matrix3x3Template operator/(RealType s) const { return Matrix3x3Template(c0 / s, c1 / s, c2 / s); }
        friend inline Matrix3x3Template operator*(RealType s, const Matrix3x3Template &mat) { return Matrix3x3Template(s * mat.c0, s * mat.c1, s * mat.c2); }
        
        Matrix3x3Template &operator+=(const Matrix3x3Template &mat) { c0 += mat.c0; c1 += mat.c1; c2 += mat.c2; return *this; }
        Matrix3x3Template &operator-=(const Matrix3x3Template &mat) { c0 -= mat.c0; c1 -= mat.c1; c2 -= mat.c2; return *this; }
        Matrix3x3Template &operator*=(const Matrix3x3Template &mat) {
            const Vector3Template<RealType> r[] = {row(0), row(1), row(2)};
            c0 = Vector3Template<RealType>(dot(r[0], mat.c0), dot(r[1], mat.c0), dot(r[2], mat.c0));
            c1 = Vector3Template<RealType>(dot(r[0], mat.c1), dot(r[1], mat.c1), dot(r[2], mat.c1));
            c2 = Vector3Template<RealType>(dot(r[0], mat.c2), dot(r[1], mat.c2), dot(r[2], mat.c2));
            return *this;
        }
        Matrix3x3Template &operator*=(RealType s) { c0 *= s; c1 *= s; c2 *= s; return *this; }
        Matrix3x3Template &operator/=(RealType s) { c0 /= s; c1 /= s; c2 /= s; return *this; }
        
        bool operator==(const Matrix3x3Template &m) const { return c0 == m.c0 && c1 == m.c1 && c2 == m.c2; }
        bool operator!=(const Matrix3x3Template &m) const { return c0 == m.c0 || c1 != m.c1 || c2 != m.c2; }
        
        Vector3Template<RealType> &operator[](unsigned int c) {
            SLRAssert(c < 3, "\"c\" is out of range [0, 2].");
            return *(&c0 + c);
        }
        
        Vector3Template<RealType> operator[](unsigned int c) const {
            SLRAssert(c < 3, "\"c\" is out of range [0, 2].");
            return *(&c0 + c);
        }
        
        const Vector3Template<RealType> &column(unsigned int c) const {
            SLRAssert(c < 3, "\"c\" is out of range [0, 2].");
            return *(&c0 + c);
        }
        Vector3Template<RealType> row(unsigned int r) const {
            SLRAssert(r < 3, "\"r\" is out of range [0, 2].");
            switch (r) {
                case 0:
                    return Vector3Template<RealType>(m00, m01, m02);
                case 1:
                    return Vector3Template<RealType>(m10, m11, m12);
                case 2:
                    return Vector3Template<RealType>(m20, m21, m22);
                default:
                    return Vector3Template<RealType>::Zero;
            }
        }
        
        Matrix3x3Template &swapColumns(unsigned int ca, unsigned int cb) {
            if (ca != cb) {
                Vector3Template<RealType> temp = column(ca);
                (*this)[ca] = (*this)[cb];
                (*this)[cb] = temp;
            }
            return *this;
        }
        
        Matrix3x3Template &swapRows(unsigned int ra, unsigned int rb) {
            if (ra != rb) {
                Vector3Template<RealType> temp = row(ra);
                setRow(ra, row(rb));
                setRow(rb, temp);
            }
            return *this;
        }
        
        Matrix3x3Template &setRow(unsigned int r, const Vector3Template<RealType> &v) {
            SLRAssert(r < 3, "\"r\" is out of range [0, 2].");
            c0[r] = v[0]; c1[r] = v[1]; c2[r] = v[2];
            return *this;
        }
        Matrix3x3Template &scaleRow(unsigned int r, RealType s) {
            SLRAssert(r < 3, "\"r\" is out of range [0, 2].");
            c0[r] *= s; c1[r] *= s; c2[r] *= s;
            return *this;
        }
        Matrix3x3Template &addRow(unsigned int r, const Vector3Template<RealType> &v) {
            SLRAssert(r < 3, "\"r\" is out of range [0, 2].");
            c0[r] += v[0]; c1[r] += v[1]; c2[r] += v[2];
            return *this;
        }
        
        RealType determinant() const {
            return (c0[0] * (c1[1] * c2[2] - c2[1] * c1[2]) -
                    c1[0] * (c0[1] * c2[2] - c2[1] * c0[2]) +
                    c2[0] * (c0[1] * c1[2] - c1[1] * c0[2]));
        }
        
        Matrix3x3Template& transpose() {
            std::swap(m10, m01); std::swap(m20, m02);
            std::swap(m21, m12);
            return *this;
        }
        
        Matrix3x3Template &invert() {
            SLRAssert_NotImplemented();
            return *this;
        }
        
        bool isIdentity() const {
            typedef Vector3Template<RealType> V3;
            return c0 == V3::Ex && c1 == V3::Ey && c2 == V3::Ez;
        }
        bool hasNaN() const { return c0.hasNaN() || c1.hasNaN() || c2.hasNaN(); }
        bool hasInf() const { return c0.hasInf() || c1.hasInf() || c2.hasInf(); }
        
        static const Matrix3x3Template Identity;
        static const Matrix3x3Template Zero;
        static const Matrix3x3Template One;
        static const Matrix3x3Template Inf;
        static const Matrix3x3Template NaN;
    };
    template <typename RealType>
    const Matrix3x3Template<RealType> Matrix3x3Template<RealType>::Identity = Matrix3x3Template<RealType>(Vector3Template<RealType>(1, 0, 0),
                                                                                                          Vector3Template<RealType>(0, 1, 0),
                                                                                                          Vector3Template<RealType>(0, 0, 1));
    template <typename RealType>
    const Matrix3x3Template<RealType> Matrix3x3Template<RealType>::Zero = Matrix3x3Template<RealType>(Vector3Template<RealType>(0),
                                                                                                      Vector3Template<RealType>(0),
                                                                                                      Vector3Template<RealType>(0));
    template <typename RealType>
    const Matrix3x3Template<RealType> Matrix3x3Template<RealType>::One = Matrix3x3Template<RealType>(Vector3Template<RealType>(1),
                                                                                                     Vector3Template<RealType>(1),
                                                                                                     Vector3Template<RealType>(1));
    template <typename RealType>
    const Matrix3x3Template<RealType> Matrix3x3Template<RealType>::Inf = Matrix3x3Template<RealType>(Vector3Template<RealType>(std::numeric_limits<RealType>::infinity()),
                                                                                                     Vector3Template<RealType>(std::numeric_limits<RealType>::infinity()),
                                                                                                     Vector3Template<RealType>(std::numeric_limits<RealType>::infinity()));
    template <typename RealType>
    const Matrix3x3Template<RealType> Matrix3x3Template<RealType>::NaN = Matrix3x3Template<RealType>(Vector3Template<RealType>(std::numeric_limits<RealType>::quiet_NaN()),
                                                                                                     Vector3Template<RealType>(std::numeric_limits<RealType>::quiet_NaN()),
                                                                                                     Vector3Template<RealType>(std::numeric_limits<RealType>::quiet_NaN()));
    
    
    template <typename RealType>
    SLR_API Matrix3x3Template<RealType> transpose(const Matrix3x3Template<RealType> &m);
    
    template <typename RealType>
    SLR_API Matrix3x3Template<RealType> invert(const Matrix3x3Template<RealType> &m);
    
    template <typename RealType>
    inline Matrix3x3Template<RealType> scale3x3(const Vector3Template<RealType> &s) {
        return Matrix3x3Template<RealType>(s.x * Vector3Template<RealType>::Ex,
                                           s.y * Vector3Template<RealType>::Ey,
                                           s.z * Vector3Template<RealType>::Ez);
    }
    
    template <typename RealType>
    inline Matrix3x3Template<RealType> scale3x3(RealType sx, RealType sy, RealType sz) {
        return scale3x3(Vector3Template<RealType>(sx, sy, sz));
    }
    
    template <typename RealType>
    inline Matrix3x3Template<RealType> scale3x3(RealType s) {
        return scale3x3(Vector3Template<RealType>(s, s, s));
    }
    
    template <typename RealType>
    SLR_API Matrix3x3Template<RealType> rotate3x3(RealType angle, const Vector3Template<RealType> &axis);
    
    template <typename RealType>
    inline Matrix3x3Template<RealType> rotate3x3(RealType angle, RealType ax, RealType ay, RealType az) {
        return rotate3x3(angle, Vector3Template<RealType>(ax, ay, az));
    }
    
    template <typename RealType>
    inline Matrix3x3Template<RealType> rotateX3x3(RealType angle) { return rotate3x3(angle, Vector3Template<RealType>::Ex); }
    template <typename RealType>
    inline Matrix3x3Template<RealType> rotateY3x3(RealType angle) { return rotate3x3(angle, Vector3Template<RealType>::Ey); }
    template <typename RealType>
    inline Matrix3x3Template<RealType> rotateZ3x3(RealType angle) { return rotate3x3(angle, Vector3Template<RealType>::Ez); }
}

#endif /* Matrix3x3_h */
