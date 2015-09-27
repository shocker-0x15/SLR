//
//  RGBTypes.h
//
//  Created by 渡部 心 on 2015/05/04.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__RGBTypes__
#define __SLR__RGBTypes__

namespace SLR {
    // just for compatibility with spectral representation build.
    template <typename RealType>
    struct RGBSamplesTemplate {
        enum Flag : uint16_t {
            LambdaSelected = 0x01,
        };
        uint16_t selectedLambda;
        uint16_t flags;
        static const uint32_t NumComponents;
        
        RGBSamplesTemplate() { };
        
        bool lambdaSelected() const {
            return (flags & LambdaSelected) != 0;
        };
        
        static RGBSamplesTemplate createWithEqualOffsets(RealType offset, RealType uLambda, RealType* PDF) {
            SLRAssert(offset >= 0 && offset < 1, "\"offset\" must be in range [0, 1).");
            SLRAssert(uLambda >= 0 && uLambda < 1, "\"uLambda\" must be in range [0, 1).");
            RGBSamplesTemplate ret;
            ret.selectedLambda = std::min(uint16_t(3 * uLambda), uint16_t(2));
            ret.flags = 0;
            *PDF = 1;
            return ret;
        };
    };
    template <typename RealType> const uint32_t RGBSamplesTemplate<RealType>::NumComponents = 3;
    
    template <typename RealType>
    struct RGBTemplate {
        RealType r, g, b;
        static const uint32_t NumComponents;
        
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
        
        RealType luminance(RGBColorSpace space = RGBColorSpace::sRGB) const {
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
        
        std::string toString() const {
            char str[256];
            sprintf(str, "(%g, %g, %g)\n", r, g, b);
            return str;
        };
        
        //------------------------------------------------
        // Methods for compatibility with ContinuousSpectrumTemplate
        const RGBTemplate &evaluate(const RGBSamplesTemplate<RealType> &wls) const {
            return *this;
        };
        //------------------------------------------------
        // Methods for compatibility with DiscretizedSpectrumTemplate
        bool hasMinus() const {
            return r < 0 || g < 0 || b < 0;
        };
        static void init() { };
        //------------------------------------------------
        
        static const RGBTemplate Zero;
        static const RGBTemplate One;
        static const RGBTemplate Inf;
        static const RGBTemplate NaN;
    };
    template <typename RealType> const uint32_t RGBTemplate<RealType>::NumComponents = 3;
    
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
        CompensatedSum<RGBTemplate<RealType>> value;
        
        RGBStorageTemplate(const RGBTemplate<RealType> &v = RGBTemplate<RealType>::Zero) :
        value(v) {};
        
        RGBStorageTemplate &add(const RGBSamplesTemplate<RealType> &wls, const RGBTemplate<RealType> &val) {
            value += val;
            return *this;
        };
    };    
}

#endif
