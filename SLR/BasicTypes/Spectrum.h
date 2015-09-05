//
//  Spectrum.h
//
//  Created by 渡部 心 on 2015/05/04.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__Spectrum__
#define __SLR__Spectrum__

#include "../defines.h"
#include "../references.h"
#include <limits>

enum class RGBColorSpace {
    sRGB = 0,
};

template <typename RealType>
struct SpectrumTemplate {
    RealType r, g, b;
    
    SpectrumTemplate(RealType v = 0.0f) : r(v), g(v), b(v) { };
    constexpr SpectrumTemplate(RealType rr, RealType gg, RealType bb) : r(rr), g(gg), b(bb) { };
    
    SpectrumTemplate operator+() const { return *this; };
    SpectrumTemplate operator-() const { return SpectrumTemplate(-r, -g, -b); };
    
    SpectrumTemplate operator+(const SpectrumTemplate &c) const { return SpectrumTemplate(r + c.r, g + c.g, b + c.b); };
    SpectrumTemplate operator-(const SpectrumTemplate &c) const { return SpectrumTemplate(r - c.r, g - c.g, b - c.b); };
    SpectrumTemplate operator*(const SpectrumTemplate &c) const { return SpectrumTemplate(r * c.r, g * c.g, b * c.b); };
    SpectrumTemplate operator/(const SpectrumTemplate &c) const { return SpectrumTemplate(r / c.r, g / c.g, b / c.b); };
    SpectrumTemplate operator*(RealType s) const { return SpectrumTemplate(r * s, g * s, b * s); };
    SpectrumTemplate operator/(RealType s) const { RealType rc = 1.0f / s; return SpectrumTemplate(r * rc, g * rc, b * rc); };
    friend inline SpectrumTemplate operator*(RealType s, const SpectrumTemplate &c) { return SpectrumTemplate(s * c.r, s * c.g, s * c.b); };
    
    SpectrumTemplate &operator+=(const SpectrumTemplate &c) { r += c.r; g += c.g; b += c.b; return *this; };
    SpectrumTemplate &operator-=(const SpectrumTemplate &c) { r -= c.r; g -= c.g; b -= c.b; return *this; };
    SpectrumTemplate &operator*=(const SpectrumTemplate &c) { r *= c.r; g *= c.g; b *= c.b; return *this; };
    SpectrumTemplate &operator/=(const SpectrumTemplate &c) { r /= c.r; g /= c.g; b /= c.b; return *this; };
    SpectrumTemplate &operator*=(RealType s) { r *= s; g *= s; b *= s; return *this; };
    SpectrumTemplate &operator/=(RealType s) { RealType rc = 1.0f / s; r *= rc; g *= rc; b *= rc; return *this; };
    
    bool operator==(const SpectrumTemplate &c) const { return r == c.r && g == c.g && b == c.b; };
    bool operator!=(const SpectrumTemplate &c) const { return r != c.r || g != c.g || b != c.b; };
    
    float &operator[](unsigned int index) {
        SLRAssert(index < 3, "\"index\" is out of range [0, 2].");
        return *(&r + index);
    };
    float operator[](unsigned int index) const {
        SLRAssert(index < 3, "\"index\" is out of range [0, 2].");
        return *(&r + index);
    };
    
    float luminance(RGBColorSpace space = RGBColorSpace::sRGB) const {
        switch(space) {
            case RGBColorSpace::sRGB:
                return 0.222485 * r + 0.716905 * g + 0.060610 * b;
        }
        return 0.0f;
    };
    
    RealType maxValue() const { return std::fmax(r, std::max(g, b)); };
    RealType minValue() const { return std::fmin(r, std::min(g, b)); };
    bool hasNonZero() const { return r != 0.0f || g != 0.0f || b != 0.0f; };
    bool hasNaN() const { using std::isnan; return isnan(r) || isnan(g) || isnan(b); };
    bool hasInf() const { using std::isinf; return isinf(r) || isinf(g) || isinf(b); };
    
    // TODO: consider which RGB color space should be used for rendering calculation and converting.
    void getRGB(float RGB[3], RGBColorSpace space = RGBColorSpace::sRGB) const {
        RGB[0] = r;
        RGB[1] = g;
        RGB[2] = b;
    };
    
    void print() const { printf("(%f, %f, %f)\n", r, g, b); };
    
    static const SpectrumTemplate Zero;
    static const SpectrumTemplate One;
    static const SpectrumTemplate Inf;
    static const SpectrumTemplate NaN;
};

template <typename RealType>
const SpectrumTemplate<RealType> SpectrumTemplate<RealType>::Zero = SpectrumTemplate<RealType>(0);
template <typename RealType>
const SpectrumTemplate<RealType> SpectrumTemplate<RealType>::One = SpectrumTemplate<RealType>(1);
template <typename RealType>
const SpectrumTemplate<RealType> SpectrumTemplate<RealType>::Inf = SpectrumTemplate<RealType>(std::numeric_limits<RealType>::infinity());
template <typename RealType>
const SpectrumTemplate<RealType> SpectrumTemplate<RealType>::NaN = SpectrumTemplate<RealType>(std::numeric_limits<RealType>::quiet_NaN());


template <typename RealType>
SpectrumTemplate<RealType> pow(const SpectrumTemplate<RealType> s, RealType p) {
    return SpectrumTemplate<RealType>(std::pow(s.r, p), std::pow(s.g, p), std::pow(s.b, p));
}

template <typename RealType>
SpectrumTemplate<RealType> exp(const SpectrumTemplate<RealType> s) {
    return SpectrumTemplate<RealType>(std::exp(s.r), std::exp(s.g), std::exp(s.b));
}

template <typename RealType>
SpectrumTemplate<RealType> inverseGammaCorrection(const SpectrumTemplate<RealType> s, RealType gamma = 2.2) {
    return SpectrumTemplate<RealType>(std::pow(s.r, gamma), std::pow(s.g, gamma), std::pow(s.b, gamma));
}


typedef SpectrumTemplate<float> Spectrum;

#endif
