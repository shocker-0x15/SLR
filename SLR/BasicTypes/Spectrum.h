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
#include "CompensatedSum.h"

enum class RGBColorSpace {
    sRGB = 0,
};

template <typename RealType>
struct RGBTemplate {
    RealType r, g, b;
    
    RGBTemplate(RealType v = 0.0f) : r(v), g(v), b(v) { };
    constexpr RGBTemplate(RealType rr, RealType gg, RealType bb) : r(rr), g(gg), b(bb) { };
    
    RGBTemplate operator+() const { return *this; };
    RGBTemplate operator-() const { return RGBTemplate(-r, -g, -b); };
    
    RGBTemplate operator+(const RGBTemplate &c) const { return RGBTemplate(r + c.r, g + c.g, b + c.b); };
    RGBTemplate operator-(const RGBTemplate &c) const { return RGBTemplate(r - c.r, g - c.g, b - c.b); };
    RGBTemplate operator*(const RGBTemplate &c) const { return RGBTemplate(r * c.r, g * c.g, b * c.b); };
    RGBTemplate operator/(const RGBTemplate &c) const { return RGBTemplate(r / c.r, g / c.g, b / c.b); };
    RGBTemplate operator*(RealType s) const { return RGBTemplate(r * s, g * s, b * s); };
    RGBTemplate operator/(RealType s) const { RealType rc = 1.0f / s; return RGBTemplate(r * rc, g * rc, b * rc); };
    friend inline RGBTemplate operator*(RealType s, const RGBTemplate &c) { return RGBTemplate(s * c.r, s * c.g, s * c.b); };
    
    RGBTemplate &operator+=(const RGBTemplate &c) { r += c.r; g += c.g; b += c.b; return *this; };
    RGBTemplate &operator-=(const RGBTemplate &c) { r -= c.r; g -= c.g; b -= c.b; return *this; };
    RGBTemplate &operator*=(const RGBTemplate &c) { r *= c.r; g *= c.g; b *= c.b; return *this; };
    RGBTemplate &operator/=(const RGBTemplate &c) { r /= c.r; g /= c.g; b /= c.b; return *this; };
    RGBTemplate &operator*=(RealType s) { r *= s; g *= s; b *= s; return *this; };
    RGBTemplate &operator/=(RealType s) { RealType rc = 1.0f / s; r *= rc; g *= rc; b *= rc; return *this; };
    
    bool operator==(const RGBTemplate &c) const { return r == c.r && g == c.g && b == c.b; };
    bool operator!=(const RGBTemplate &c) const { return r != c.r || g != c.g || b != c.b; };
    
    float &operator[](unsigned int index) {
        SLRAssert(index < 3, "\"index\" is out of range [0, 2].");
        return *(&r + index);
    };
    float operator[](unsigned int index) const {
        SLRAssert(index < 3, "\"index\" is out of range [0, 2].");
        return *(&r + index);
    };
    
    RealType maxValue() const { return std::fmax(r, std::max(g, b)); };
    RealType minValue() const { return std::fmin(r, std::min(g, b)); };
    bool hasNonZero() const { return r != 0.0f || g != 0.0f || b != 0.0f; };
    bool hasNaN() const { using std::isnan; return isnan(r) || isnan(g) || isnan(b); };
    bool hasInf() const { using std::isinf; return isinf(r) || isinf(g) || isinf(b); };
    
    std::string toString() const {
        char str[256];
        sprintf(str, "(%g, %g, %g)\n", r, g, b);
        return str;
    };
    
    static const RGBTemplate Zero;
    static const RGBTemplate One;
    static const RGBTemplate Inf;
    static const RGBTemplate NaN;
};

template <typename RealType>
const RGBTemplate<RealType> RGBTemplate<RealType>::Zero = RGBTemplate<RealType>(0);
template <typename RealType>
const RGBTemplate<RealType> RGBTemplate<RealType>::One = RGBTemplate<RealType>(1);
template <typename RealType>
const RGBTemplate<RealType> RGBTemplate<RealType>::Inf = RGBTemplate<RealType>(std::numeric_limits<RealType>::infinity());
template <typename RealType>
const RGBTemplate<RealType> RGBTemplate<RealType>::NaN = RGBTemplate<RealType>(std::numeric_limits<RealType>::quiet_NaN());


template <typename RealType>
RGBTemplate<RealType> pow(const RGBTemplate<RealType> s, RealType p) {
    return RGBTemplate<RealType>(std::pow(s.r, p), std::pow(s.g, p), std::pow(s.b, p));
}

template <typename RealType>
RGBTemplate<RealType> exp(const RGBTemplate<RealType> s) {
    return RGBTemplate<RealType>(std::exp(s.r), std::exp(s.g), std::exp(s.b));
}

template <typename RealType>
RGBTemplate<RealType> inverseGammaCorrection(const RGBTemplate<RealType> s, RealType gamma = 2.2) {
    return RGBTemplate<RealType>(std::pow(s.r, gamma), std::pow(s.g, gamma), std::pow(s.b, gamma));
}


template <typename RealType>
struct RGBStorageTemplate {
    RealType r, g, b;
public:
    RGBStorageTemplate(RealType v = 0) : r(v), g(v), b(v) { };
    constexpr RGBStorageTemplate(RealType rr, RealType gg, RealType bb) : r(rr), g(gg), b(bb) { };
    
    RGBStorageTemplate operator+() const { return *this; };
    RGBStorageTemplate operator-() const { return RGBStorageTemplate(-r, -g, -b); };
    
    RGBStorageTemplate operator+(const RGBStorageTemplate &c) const { return RGBStorageTemplate(r + c.r, g + c.g, b + c.b); };
    RGBStorageTemplate operator-(const RGBStorageTemplate &c) const { return RGBStorageTemplate(r - c.r, g - c.g, b - c.b); };
    RGBStorageTemplate operator*(const RGBStorageTemplate &c) const { return RGBStorageTemplate(r * c.r, g * c.g, b * c.b); };
    RGBStorageTemplate operator/(const RGBStorageTemplate &c) const { return RGBStorageTemplate(r / c.r, g / c.g, b / c.b); };
    RGBStorageTemplate operator*(RealType s) const { return RGBStorageTemplate(r * s, g * s, b * s); };
    RGBStorageTemplate operator/(RealType s) const { RealType rc = 1.0f / s; return RGBStorageTemplate(r * rc, g * rc, b * rc); };
    friend inline RGBStorageTemplate operator*(RealType s, const RGBStorageTemplate &c) { return RGBStorageTemplate(s * c.r, s * c.g, s * c.b); };
    
    RGBStorageTemplate &operator+=(const RGBStorageTemplate &c) { r += c.r; g += c.g; b += c.b; return *this; };
    RGBStorageTemplate &operator-=(const RGBStorageTemplate &c) { r -= c.r; g -= c.g; b -= c.b; return *this; };
    RGBStorageTemplate &operator*=(const RGBStorageTemplate &c) { r *= c.r; g *= c.g; b *= c.b; return *this; };
    RGBStorageTemplate &operator/=(const RGBStorageTemplate &c) { r /= c.r; g /= c.g; b /= c.b; return *this; };
    RGBStorageTemplate &operator*=(RealType s) { r *= s; g *= s; b *= s; return *this; };
    RGBStorageTemplate &operator/=(RealType s) { RealType rc = 1.0f / s; r *= rc; g *= rc; b *= rc; return *this; };
    
    bool operator==(const RGBStorageTemplate &c) const { return r == c.r && g == c.g && b == c.b; };
    bool operator!=(const RGBStorageTemplate &c) const { return r != c.r || g != c.g || b != c.b; };
    
    RealType &operator[](unsigned int index) {
        SLRAssert(index < 3, "\"index\" is out of range [0, 2].");
        return *(&r + index);
    };
    RealType operator[](unsigned int index) const {
        SLRAssert(index < 3, "\"index\" is out of range [0, 2].");
        return *(&r + index);
    };
    
    RealType maxValue() const { return std::fmax(r, std::max(g, b)); };
    RealType minValue() const { return std::fmin(r, std::min(g, b)); };
    bool hasNonZero() const { return r != 0.0f || g != 0.0f || b != 0.0f; };
    bool hasNaN() const { using std::isnan; return isnan(r) || isnan(g) || isnan(b); };
    bool hasInf() const { using std::isinf; return isinf(r) || isinf(g) || isinf(b); };
    
    std::string toString() const {
        char str[256];
        sprintf(str, "(%g, %g, %g)\n", r, g, b);
        return str;
    };
    
    RGBStorageTemplate &add(const RGBTemplate<RealType> &value) {
        r += value.r;
        g += value.g;
        b += value.b;
        return *this;
    };
    
    float luminance(RGBColorSpace space = RGBColorSpace::sRGB) const {
        switch(space) {
            case RGBColorSpace::sRGB:
                return 0.222485 * r + 0.716905 * g + 0.060610 * b;
        }
        return 0.0f;
    };
    
    // TODO: consider which RGB color space should be used for rendering calculation and converting.
    void getRGB(float RGB[3], RGBColorSpace space = RGBColorSpace::sRGB) const {
        RGB[0] = r;
        RGB[1] = g;
        RGB[2] = b;
    };
};


template <typename RealType, uint32_t N>
struct WavelengthSamplesTemplate {
    RealType lambda[N];
};

template <typename RealType, uint32_t N>
struct SampledSpectrumTemplate {
    RealType values[N];
    static const uint32_t numComponents;
    
    SampledSpectrumTemplate(RealType v = 0.0f) { for (int i = 0; i < N; ++i) values[i] = v; };
    SampledSpectrumTemplate(const RealType* vals) { for (int i = 0; i < N; ++i) values[i] = vals[i]; };
    
    SampledSpectrumTemplate operator+() const { return *this; };
    SampledSpectrumTemplate operator-() const {
        RealType vals[N];
        for (int i = 0; i < N; ++i)
            vals[i] = -values[i];
        return SampledSpectrumTemplate(vals);
    };
    
    SampledSpectrumTemplate operator+(const SampledSpectrumTemplate &c) const {
        RealType vals[N];
        for (int i = 0; i < N; ++i)
            vals[i] = values[i] + c.values[i];
        return SampledSpectrumTemplate(vals);
    };
    SampledSpectrumTemplate operator-(const SampledSpectrumTemplate &c) const {
        RealType vals[N];
        for (int i = 0; i < N; ++i)
            vals[i] = values[i] - c.values[i];
        return SampledSpectrumTemplate(vals);
    };
    SampledSpectrumTemplate operator*(const SampledSpectrumTemplate &c) const {
        RealType vals[N];
        for (int i = 0; i < N; ++i)
            vals[i] = values[i] * c.values[i];
        return SampledSpectrumTemplate(vals);
    };
    SampledSpectrumTemplate operator/(const SampledSpectrumTemplate &c) const {
        RealType vals[N];
        for (int i = 0; i < N; ++i)
            vals[i] = values[i] / c.values[i];
        return SampledSpectrumTemplate(vals);
    };
    SampledSpectrumTemplate operator*(RealType s) const {
        RealType vals[N];
        for (int i = 0; i < N; ++i)
            vals[i] = values[i] * s;
        return SampledSpectrumTemplate(vals);
    };
    SampledSpectrumTemplate operator/(RealType s) const {
        RealType vals[N];
        RealType r = 1 / s;
        for (int i = 0; i < N; ++i)
            vals[i] = values[i] * r;
        return SampledSpectrumTemplate(vals);
    };
    friend inline SampledSpectrumTemplate operator*(RealType s, const SampledSpectrumTemplate &c) {
        RealType vals[N];
        for (int i = 0; i < N; ++i)
            vals[i] = c.values[i] * s;
        return SampledSpectrumTemplate(vals);
    };
    
    SampledSpectrumTemplate &operator+=(const SampledSpectrumTemplate &c) {
        for (int i = 0; i < N; ++i)
            values[i] += c.values[i];
        return *this;
    };
    SampledSpectrumTemplate &operator-=(const SampledSpectrumTemplate &c) {
        for (int i = 0; i < N; ++i)
            values[i] -= c.values[i];
        return *this;
    };
    SampledSpectrumTemplate &operator*=(const SampledSpectrumTemplate &c) {
        for (int i = 0; i < N; ++i)
            values[i] *= c.values[i];
        return *this;
    };
    SampledSpectrumTemplate &operator/=(const SampledSpectrumTemplate &c) {
        for (int i = 0; i < N; ++i)
            values[i] /= c.values[i];
        return *this;
    };
    SampledSpectrumTemplate &operator*=(RealType s) {
        for (int i = 0; i < N; ++i)
            values[i] *= s;
        return *this;
    };
    SampledSpectrumTemplate &operator/=(RealType s) {
        float r = 1 / s;
        for (int i = 0; i < N; ++i)
            values[i] *= r;
        return *this;
    };
    
    bool operator==(const SampledSpectrumTemplate &c) const {
        for (int i = 0; i < N; ++i)
            if (values[i] != c.values[i])
                return false;
        return true;
    };
    bool operator!=(const SampledSpectrumTemplate &c) const {
        for (int i = 0; i < N; ++i)
            if (values[i] != c.values[i])
                return true;
        return false;
    };
    
    RealType &operator[](unsigned int index) {
        SLRAssert(index < N, "\"index\" is out of range [0, %u].", N);
        return values[index];
    };
    RealType operator[](unsigned int index) const {
        SLRAssert(index < N, "\"index\" is out of range [0, %u].", N);
        return values[index];
    };
    
    RealType maxValue() const {
        RealType maxVal = values[0];
        for (int i = 1; i < N; ++i)
            maxVal = std::fmax(values[i], maxVal);
        return maxVal;
    };
    RealType minValue() const {
        RealType minVal = values[0];
        for (int i = 1; i < N; ++i)
            minVal = std::fmin(values[i], minVal);
        return minVal;
    };
    bool hasNonZero() const {
        for (int i = 0; i < N; ++i)
            if (values[i] != 0)
                return true;
        return false;
    };
    bool hasNaN() const {
        for (int i = 0; i < N; ++i)
            if (std::isnan(values[i]))
                return true;
        return false;
    };
    bool hasInf() const {
        for (int i = 0; i < N; ++i)
            if (std::isinf(values[i]))
                return true;
        return false;
    };
    
    std::string toString() const {
        std::string ret = "(";
        char str[256];
        for (int i = 0; i < N - 1; ++i) {
            sprintf(str, "%g, ", values[i]);
            ret += str;
        }
        sprintf(str, "%g)", values[N - 1]);
        ret += str;
        return str;
    };
    
    static const SampledSpectrumTemplate Zero;
    static const SampledSpectrumTemplate One;
    static const SampledSpectrumTemplate Inf;
    static const SampledSpectrumTemplate NaN;
};
template <typename RealType, uint32_t N> const uint32_t SampledSpectrumTemplate<RealType, N>::numComponents = N;

template <typename RealType, uint32_t N>
const SampledSpectrumTemplate<RealType, N> SampledSpectrumTemplate<RealType, N>::Zero = SampledSpectrumTemplate<RealType, N>(0);
template <typename RealType, uint32_t N>
const SampledSpectrumTemplate<RealType, N> SampledSpectrumTemplate<RealType, N>::One = SampledSpectrumTemplate<RealType, N>(1);
template <typename RealType, uint32_t N>
const SampledSpectrumTemplate<RealType, N> SampledSpectrumTemplate<RealType, N>::Inf = SampledSpectrumTemplate<RealType, N>(std::numeric_limits<RealType>::infinity());
template <typename RealType, uint32_t N>
const SampledSpectrumTemplate<RealType, N> SampledSpectrumTemplate<RealType, N>::NaN = SampledSpectrumTemplate<RealType, N>(std::numeric_limits<RealType>::quiet_NaN());


template <typename RealType, uint32_t numStrata>
struct SpectrumStorageTemplate {
    RealType value[numStrata];

    SpectrumStorageTemplate(RealType v = 0) {
        for (int i = 0; i < numStrata; ++i)
            value[i] = v;
    };
    
    template <uint32_t N>
    SpectrumStorageTemplate &add(const WavelengthSamplesTemplate<RealType, N> &wls, const SampledSpectrumTemplate<RealType, N> &value) {
        return *this;
    };
    
    float luminance(RGBColorSpace space = RGBColorSpace::sRGB) const {
        return 0.0f;
    };
    
    void getRGB(float RGB[3], RGBColorSpace space = RGBColorSpace::sRGB) const {
    };
};

#endif
