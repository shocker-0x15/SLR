//
//  TexCoord2.h
//
//  Created by 渡部 心 on 2015/06/11.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_TexCoord2_h__
#define __SLR_TexCoord2_h__

#include "../defines.h"
#include "../references.h"

namespace SLR {
    template <typename RealType>
    struct TexCoord2Template {
        RealType u, v;
        
        TexCoord2Template(RealType val = 0.0f) : u(val), v(val) { };
        CONSTEXPR_CONSTRUCTOR TexCoord2Template(RealType uu, RealType vv) : u(uu), v(vv) { };
        
        TexCoord2Template operator+() const { return *this; };
        TexCoord2Template operator-() const { return TexCoord2Template(-u, -v); };
        
        TexCoord2Template operator+(const TexCoord2Template<RealType> &t) const { return TexCoord2Template(u + t.u, v + t.v); };
        TexCoord2Template operator-(const TexCoord2Template<RealType> &t) const { return TexCoord2Template(u - t.u, v - t.v); };
        TexCoord2Template operator*(RealType s) const { return TexCoord2Template(u * s, v * s); };
        TexCoord2Template operator/(RealType s) const { RealType r = 1.0f / s; return TexCoord2Template(u * r, v * r); };
        friend inline TexCoord2Template operator*(RealType s, const TexCoord2Template &t) { return TexCoord2Template(s * t.u, s * t.v); };
        
        TexCoord2Template &operator+=(const TexCoord2Template<RealType> &t) { u += t.u; v += t.v; return *this; };
        TexCoord2Template &operator-=(const TexCoord2Template<RealType> &t) { u -= t.u; v -= t.v; return *this; };
        TexCoord2Template &operator*=(RealType s) { u *= s; v *= s; return *this; };
        TexCoord2Template &operator/=(RealType s) { RealType r = 1.0f / s; u *= r; v *= r; return *this; };
        
        bool operator==(const TexCoord2Template &t) const { return u == t.u && v == t.v; };
        bool operator!=(const TexCoord2Template &t) const { return u != t.u || v != t.v; };
        
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
        
        static const TexCoord2Template Zero;
        static const TexCoord2Template One;
        static const TexCoord2Template Inf;
        static const TexCoord2Template NaN;
    };
    
    template <typename RealType>
    const TexCoord2Template<RealType> TexCoord2Template<RealType>::Zero = TexCoord2Template<RealType>(0);
    template <typename RealType>
    const TexCoord2Template<RealType> TexCoord2Template<RealType>::One = TexCoord2Template<RealType>(1);
    template <typename RealType>
    const TexCoord2Template<RealType> TexCoord2Template<RealType>::Inf = TexCoord2Template<RealType>(std::numeric_limits<RealType>::infinity());
    template <typename RealType>
    const TexCoord2Template<RealType> TexCoord2Template<RealType>::NaN = TexCoord2Template<RealType>(std::numeric_limits<RealType>::quiet_NaN());    
}

#endif
