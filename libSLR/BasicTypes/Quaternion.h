//
//  Quaternion.h
//
//  Created by 渡部 心 on 2015/07/20.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_Quaternion__
#define __SLR_Quaternion__

#include "../defines.h"
#include "Vector3.h"
#include "Matrix4x4.h"

namespace SLR {
    template <typename RealType>
    struct SLR_API QuaternionTemplate {
        union {
            Vector3Template<RealType> v;
            struct { RealType x, y, z; };
        };
        RealType w;
        
        QuaternionTemplate() : w(1.0f) { };
        CONSTEXPR_CONSTRUCTOR QuaternionTemplate(RealType xx, RealType yy, RealType zz, RealType ww) : v(xx, yy, zz), w(ww) { };
        CONSTEXPR_CONSTRUCTOR QuaternionTemplate(const Vector3Template<RealType> &vv, RealType ww) : v(vv), w(ww) { };
        QuaternionTemplate(const Matrix4x4Template<RealType> &m) {
            RealType trace = m[0][0] + m[1][1] + m[2][2];
            if (trace > 0.0f) {
                RealType s = std::sqrt(trace + 1.0f);
                v = (0.5f / s) * Vector3Template<RealType>(m[1][2] - m[2][1], m[2][0] - m[0][2], m[0][1] - m[1][0]);
                w = s / 2.0f;
            }
            else {
                const int nxt[3] = {1, 2, 0};
                RealType q[3];
                int i = 0;
                if (m[1][1] > m[0][0])
                    i = 1;
                if (m[2][2] > m[i][i])
                    i = 2;
                int j = nxt[i];
                int k = nxt[j];
                RealType s = std::sqrt((m[i][i] - (m[j][j] + m[k][k])) + 1.0f);
                q[i] = s * 0.5f;
                if (s != 0.0f)
                    s = 0.5f / s;
                w = (m[j][k] - m[k][j]) * s;
                q[j] = (m[i][j] + m[j][i]) * s;
                q[k] = (m[i][k] + m[k][i]) * s;
                v = Vector3Template<RealType>(q[0], q[1], q[2]);
            }
        };
        
        QuaternionTemplate operator+() const { return *this; };
        QuaternionTemplate operator-() const { return QuaternionTemplate(-v, -w); };
        
        QuaternionTemplate operator+(const QuaternionTemplate &q) const { return QuaternionTemplate(v + q.v, w + q.w); };
        QuaternionTemplate operator-(const QuaternionTemplate &q) const { return QuaternionTemplate(v - q.v, w - q.w); };
        QuaternionTemplate operator*(const QuaternionTemplate &q) const {
            return QuaternionTemplate(cross(v, q.v) + w * q.v + q.w * v, w * q.w - dot(v, q.v));
        };
        QuaternionTemplate operator*(RealType s) const { return QuaternionTemplate(v * s, w * s); };
        QuaternionTemplate operator/(RealType s) const { RealType r = 1.0f / s; return QuaternionTemplate(v * r, w * r); };
        friend inline QuaternionTemplate operator*(RealType s, const QuaternionTemplate &q) { return QuaternionTemplate(q.v * s, q.w * s); };
        
        QuaternionTemplate &operator+=(const QuaternionTemplate &q) { v += q.v; w += q.w; return *this; };
        QuaternionTemplate &operator-=(const QuaternionTemplate &q) { v -= q.v; w -= q.w; return *this; };
        QuaternionTemplate &operator*=(RealType s) { v *= s; w *= s; return *this; };
        QuaternionTemplate &operator/=(RealType s) { RealType r = 1.0f / s; v *= r; w *= r; return *this; };
        
        Matrix4x4Template<RealType> toMatrix() const {
            RealType xx = x * x, yy = y * y, zz = z * z;
            RealType xy = x * y, yz = y * z, zx = z * x;
            RealType xw = x * w, yw = y * w, zw = z * w;
            return Matrix4x4Template<RealType>(Vector3Template<RealType>(1 - 2 * (yy + zz), 2 * (xy + zw), 2 * (zx - yw)),
                                               Vector3Template<RealType>(2 * (xy - zw), 1 - 2 * (xx + zz), 2 * (yz + xw)),
                                               Vector3Template<RealType>(2 * (zx + yw), 2 * (yz - xw), 1 - 2 * (xx + yy)));
        };
        
        bool operator==(const QuaternionTemplate &q) const { return v == q.v && w == q.w; };
        bool operator!=(const QuaternionTemplate &q) const { return v != q.v || w != q.w; };
    };
    
    
    template <typename RealType>
    inline RealType dot(const QuaternionTemplate<RealType> &q0, const QuaternionTemplate<RealType> &q1) {
        return dot(q0.v, q1.v) + q0.w * q1.w;
    }
    
    template <typename RealType>
    inline QuaternionTemplate<RealType> normalize(const QuaternionTemplate<RealType> &q) {
        return q / std::sqrt(dot(q, q));
    }
    
    template <typename RealType>
    inline QuaternionTemplate<RealType> Slerp(RealType t, const QuaternionTemplate<RealType> &q0, const QuaternionTemplate<RealType> &q1) {
        RealType cosTheta = dot(q0, q1);
        if (cosTheta > 0.9995f)
            return normalize((1 - t) * q0 + t * q1);
        else {
            RealType theta = std::acos(std::clamp(cosTheta, (RealType)-1, (RealType)1));
            RealType thetap = theta * t;
            QuaternionTemplate<RealType> qPerp = normalize(q1 - q0 * cosTheta);
            return q0 * std::cos(thetap) + qPerp * std::sin(thetap);
        }
    }
    
    template <typename RealType>
    SLR_API void decompose(const Matrix4x4Template<RealType> &mat, Vector3Template<RealType>* T, QuaternionTemplate<RealType>* R, Matrix4x4Template<RealType>* S);
}

#endif /* __SLR_Quaternion__ */
