//
//  TexCoord2D.h
//
//  Created by 渡部 心 on 2015/06/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_TexCoord2__
#define __SLR_TexCoord2__

#include "../defines.h"
#include "../declarations.h"

namespace SLR {
    template <typename RealType>
    struct SLR_API TexCoord2DTemplate {
        RealType u, v;
        
        TexCoord2DTemplate(RealType val = 0.0f) : u(val), v(val) { };
        CONSTEXPR_CONSTRUCTOR TexCoord2DTemplate(RealType uu, RealType vv) : u(uu), v(vv) { };
        
        TexCoord2DTemplate operator+() const { return *this; };
        TexCoord2DTemplate operator-() const { return TexCoord2DTemplate(-u, -v); };
        
        TexCoord2DTemplate operator+(const TexCoord2DTemplate<RealType> &t) const { return TexCoord2DTemplate(u + t.u, v + t.v); };
        TexCoord2DTemplate operator-(const TexCoord2DTemplate<RealType> &t) const { return TexCoord2DTemplate(u - t.u, v - t.v); };
        TexCoord2DTemplate operator*(RealType s) const { return TexCoord2DTemplate(u * s, v * s); };
        TexCoord2DTemplate operator/(RealType s) const { RealType r = 1.0f / s; return TexCoord2DTemplate(u * r, v * r); };
        friend inline TexCoord2DTemplate operator*(RealType s, const TexCoord2DTemplate &t) { return TexCoord2DTemplate(s * t.u, s * t.v); };
        
        TexCoord2DTemplate &operator+=(const TexCoord2DTemplate<RealType> &t) { u += t.u; v += t.v; return *this; };
        TexCoord2DTemplate &operator-=(const TexCoord2DTemplate<RealType> &t) { u -= t.u; v -= t.v; return *this; };
        TexCoord2DTemplate &operator*=(RealType s) { u *= s; v *= s; return *this; };
        TexCoord2DTemplate &operator/=(RealType s) { RealType r = 1.0f / s; u *= r; v *= r; return *this; };
        
        bool operator==(const TexCoord2DTemplate &t) const { return u == t.u && v == t.v; };
        bool operator!=(const TexCoord2DTemplate &t) const { return u != t.u || v != t.v; };
        
        RealType &operator[](unsigned int index) {
            SLRAssert(index < 2, "\"index\" is out of range [0, 1].");
            return *(&u + index);
        };
        RealType operator[](unsigned int index) const {
            SLRAssert(index < 2, "\"index\" is out of range [0, 1].");
            return *(&u + index);
        };
        
        RealType maxValue() const { return std::fmax(u, v); };
        RealType minValue() const { return std::fmin(u, v); };
        bool hasNaN() const { using std::isnan; return isnan(u) || isnan(v); };
        bool hasInf() const { using std::isinf; return isinf(u) || isinf(v); };
        
        void print() const { printf("(%g, %g)\n", u, v); };
        
        static const TexCoord2DTemplate Zero;
        static const TexCoord2DTemplate One;
        static const TexCoord2DTemplate Inf;
        static const TexCoord2DTemplate NaN;
    };  
    template <typename RealType>
    const TexCoord2DTemplate<RealType> TexCoord2DTemplate<RealType>::Zero = TexCoord2DTemplate<RealType>(0);
    template <typename RealType>
    const TexCoord2DTemplate<RealType> TexCoord2DTemplate<RealType>::One = TexCoord2DTemplate<RealType>(1);
    template <typename RealType>
    const TexCoord2DTemplate<RealType> TexCoord2DTemplate<RealType>::Inf = TexCoord2DTemplate<RealType>(std::numeric_limits<RealType>::infinity());
    template <typename RealType>
    const TexCoord2DTemplate<RealType> TexCoord2DTemplate<RealType>::NaN = TexCoord2DTemplate<RealType>(std::numeric_limits<RealType>::quiet_NaN());
}

#endif /* __SLR_TexCoord2__ */
