//
//  SpectrumTypes.h
//
//  Created by 渡部 心 on 2015/09/13.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLR_SpectrumTypes_h
#define SLR_SpectrumTypes_h

#include "../defines.h"
#include "../references.h"
#include "Spectrum.h"
#include "CompensatedSum.h"

namespace SLR {
    template <typename RealType, uint32_t N>
    struct SLR_API WavelengthSamplesTemplate {
        enum Flag : uint16_t {
            LambdaIsSelected = 0x01,
        };

        RealType lambdas[N];
        uint16_t selectedLambda;
        uint16_t flags;
        
        WavelengthSamplesTemplate() : selectedLambda(0), flags(0) {};
        WavelengthSamplesTemplate(const RealType* values) {
            for (int i = 0; i < N; ++i)
                lambdas[i] = values[i];
            selectedLambda = 0;
            flags = 0;
        }
        WavelengthSamplesTemplate(const WavelengthSamplesTemplate &wls) {
            for (int i = 0; i < N; ++i)
                lambdas[i] = wls.lambdas[i];
            selectedLambda = wls.selectedLambda;
            flags = wls.flags;
        }
        
        RealType &operator[](uint32_t index) {
            SLRAssert(index < N, "\"index\" is out of range [0, %u].", N - 1);
            return lambdas[index];
        }
        RealType operator[](uint32_t index) const {
            SLRAssert(index < N, "\"index\" is out of range [0, %u].", N - 1);
            return lambdas[index];
        }
        
        bool lambdaSelected() const {
            return (flags & LambdaIsSelected) != 0;
        }
        
        static WavelengthSamplesTemplate createWithEqualOffsets(RealType offset, RealType uLambda, RealType* PDF) {
            SLRAssert(offset >= 0 && offset < 1, "\"offset\" must be in range [0, 1).");
            SLRAssert(uLambda >= 0 && uLambda < 1, "\"uLambda\" must be in range [0, 1).");
            WavelengthSamplesTemplate wls;
            for (int i = 0; i < N; ++i)
                wls.lambdas[i] = WavelengthLowBound + (WavelengthHighBound - WavelengthLowBound) * (i + offset) / N;
            wls.selectedLambda = std::min(uint16_t(N * uLambda), uint16_t(N - 1));
            wls.flags = 0;
            *PDF = N / (WavelengthHighBound - WavelengthLowBound);
            return wls;
        }
        
        static const uint32_t NumComponents;
    };
    template <typename RealType, uint32_t N>
    const uint32_t WavelengthSamplesTemplate<RealType, N>::NumComponents = N;

    
    
    template <typename RealType, uint32_t N>
    class SLR_API ContinuousSpectrumTemplate {
    public:
        virtual ~ContinuousSpectrumTemplate() { }
        virtual RealType calcBounds() const = 0;
        virtual SampledSpectrumTemplate<RealType, N> evaluate(const WavelengthSamplesTemplate<RealType, N> &wls) const = 0;
        virtual ContinuousSpectrumTemplate* createScaled(RealType scale) const = 0;
    };
    
    
    
    template <typename RealType, uint32_t N>
    class SLR_API RegularContinuousSpectrumTemplate : public ContinuousSpectrumTemplate<RealType, N> {
        RealType minLambda, maxLambda;
        uint32_t numSamples;
        RealType* values;
        
    public:
        RegularContinuousSpectrumTemplate(RealType minWL, RealType maxWL, const RealType* vals, uint32_t numVals) : minLambda(minWL), maxLambda(maxWL), numSamples(numVals) {
            values = new RealType[numSamples];
            for (int i = 0; i < numSamples; ++i)
                values[i] = vals[i];
        }
        ~RegularContinuousSpectrumTemplate() {
            delete[] values;
        }
        
        RealType calcBounds() const override;
        
        SampledSpectrumTemplate<RealType, N> evaluate(const WavelengthSamplesTemplate<RealType, N> &wls) const override;
        
        ContinuousSpectrumTemplate<RealType, N>* createScaled(RealType scale) const override;
    };

    
    
    template <typename RealType, uint32_t N>
    class SLR_API IrregularContinuousSpectrumTemplate : public ContinuousSpectrumTemplate<RealType, N> {
        uint32_t numSamples;
        RealType* lambdas;
        RealType* values;
        
    public:
        IrregularContinuousSpectrumTemplate(const RealType* wls, const RealType* vals, uint32_t numVals) : numSamples(numVals) {
            lambdas = new RealType[numSamples];
            values = new RealType[numSamples];
            for (int i = 0; i < numSamples; ++i) {
                lambdas[i] = wls[i];
                values[i] = vals[i];
            }
        }
        ~IrregularContinuousSpectrumTemplate() {
            delete[] lambdas;
            delete[] values;
        }
        
        RealType calcBounds() const override;
        
        SampledSpectrumTemplate<RealType, N> evaluate(const WavelengthSamplesTemplate<RealType, N> &wls) const override;
        
        ContinuousSpectrumTemplate<RealType, N>* createScaled(RealType scale) const override;
    };

    
    
    // References
    // Physically Meaningful Rendering using Tristimulus Colours
    template <typename RealType, uint32_t N>
    class SLR_API UpsampledContinuousSpectrumTemplate : public ContinuousSpectrumTemplate<RealType, N> {
        uint32_t m_adjIndices;
        RealType m_s, m_t;
        RealType m_scale;
        
        void computeAdjacents(RealType u, RealType v);
        
    public:
        UpsampledContinuousSpectrumTemplate(uint32_t adjIndices, RealType s, RealType t, RealType scale) :
        m_adjIndices(adjIndices), m_s(s), m_t(t), m_scale(scale) { }
        
        UpsampledContinuousSpectrumTemplate(RealType u, RealType v, RealType scale) {
            computeAdjacents(u, v);
            m_scale = scale;
        }
        
        UpsampledContinuousSpectrumTemplate(SpectrumType spType, ColorSpace space, RealType e0, RealType e1, RealType e2);
        
        RealType calcBounds() const override;
        
        SampledSpectrumTemplate<RealType, N> evaluate(const WavelengthSamplesTemplate<RealType, N> &wls) const override;
        
        ContinuousSpectrumTemplate<RealType, N>* createScaled(RealType scale) const override {
            return new UpsampledContinuousSpectrumTemplate(m_adjIndices, m_s, m_t, this->m_scale * scale);
        }
    };

    

    template <typename RealType, uint32_t N>
    struct SLR_API SampledSpectrumTemplate {
        RealType values[N];

        SampledSpectrumTemplate(RealType v = 0.0f) { for (int i = 0; i < N; ++i) values[i] = v; }
        SampledSpectrumTemplate(const RealType* vals) { for (int i = 0; i < N; ++i) values[i] = vals[i]; }
        
        SampledSpectrumTemplate operator+() const { return *this; };
        SampledSpectrumTemplate operator-() const {
            RealType vals[N];
            for (int i = 0; i < N; ++i)
                vals[i] = -values[i];
            return SampledSpectrumTemplate(vals);
        }
        
        SampledSpectrumTemplate operator+(const SampledSpectrumTemplate &c) const {
            RealType vals[N];
            for (int i = 0; i < N; ++i)
                vals[i] = values[i] + c.values[i];
            return SampledSpectrumTemplate(vals);
        }
        SampledSpectrumTemplate operator-(const SampledSpectrumTemplate &c) const {
            RealType vals[N];
            for (int i = 0; i < N; ++i)
                vals[i] = values[i] - c.values[i];
            return SampledSpectrumTemplate(vals);
        }
        SampledSpectrumTemplate operator*(const SampledSpectrumTemplate &c) const {
            RealType vals[N];
            for (int i = 0; i < N; ++i)
                vals[i] = values[i] * c.values[i];
            return SampledSpectrumTemplate(vals);
        }
        SampledSpectrumTemplate operator/(const SampledSpectrumTemplate &c) const {
            RealType vals[N];
            for (int i = 0; i < N; ++i)
                vals[i] = values[i] / c.values[i];
            return SampledSpectrumTemplate(vals);
        }
        SampledSpectrumTemplate operator*(RealType s) const {
            RealType vals[N];
            for (int i = 0; i < N; ++i)
                vals[i] = values[i] * s;
            return SampledSpectrumTemplate(vals);
        }
        SampledSpectrumTemplate operator/(RealType s) const {
            RealType vals[N];
            RealType r = 1 / s;
            for (int i = 0; i < N; ++i)
                vals[i] = values[i] * r;
            return SampledSpectrumTemplate(vals);
        }
        friend inline SampledSpectrumTemplate operator*(RealType s, const SampledSpectrumTemplate &c) {
            RealType vals[N];
            for (int i = 0; i < N; ++i)
                vals[i] = c.values[i] * s;
            return SampledSpectrumTemplate(vals);
        }
        
        SampledSpectrumTemplate &operator+=(const SampledSpectrumTemplate &c) {
            for (int i = 0; i < N; ++i)
                values[i] += c.values[i];
            return *this;
        }
        SampledSpectrumTemplate &operator-=(const SampledSpectrumTemplate &c) {
            for (int i = 0; i < N; ++i)
                values[i] -= c.values[i];
            return *this;
        }
        SampledSpectrumTemplate &operator*=(const SampledSpectrumTemplate &c) {
            for (int i = 0; i < N; ++i)
                values[i] *= c.values[i];
            return *this;
        }
        SampledSpectrumTemplate &operator/=(const SampledSpectrumTemplate &c) {
            for (int i = 0; i < N; ++i)
                values[i] /= c.values[i];
            return *this;
        }
        SampledSpectrumTemplate &operator*=(RealType s) {
            for (int i = 0; i < N; ++i)
                values[i] *= s;
            return *this;
        }
        SampledSpectrumTemplate &operator/=(RealType s) {
            RealType r = 1 / s;
            for (int i = 0; i < N; ++i)
                values[i] *= r;
            return *this;
        }
        
        bool operator==(const SampledSpectrumTemplate &c) const {
            for (int i = 0; i < N; ++i)
                if (values[i] != c.values[i])
                    return false;
            return true;
        }
        bool operator!=(const SampledSpectrumTemplate &c) const {
            for (int i = 0; i < N; ++i)
                if (values[i] != c.values[i])
                    return true;
            return false;
        }
        
        RealType &operator[](unsigned int index) {
            SLRAssert(index < N, "\"index\" is out of range [0, %u].", N - 1);
            return values[index];
        }
        RealType operator[](unsigned int index) const {
            SLRAssert(index < N, "\"index\" is out of range [0, %u].", N - 1);
            return values[index];
        }
        
        RealType avgValue() const {
            RealType sumVal = values[0];
            for (int i = 1; i < N; ++i)
                sumVal += values[i];
            return sumVal / N;
        }
        RealType maxValue() const {
            RealType maxVal = values[0];
            for (int i = 1; i < N; ++i)
                maxVal = std::fmax(values[i], maxVal);
            return maxVal;
        }
        RealType minValue() const {
            RealType minVal = values[0];
            for (int i = 1; i < N; ++i)
                minVal = std::fmin(values[i], minVal);
            return minVal;
        }
        bool hasNonZero() const {
            for (int i = 0; i < N; ++i)
                if (values[i] != 0)
                    return true;
            return false;
        }
        bool hasNaN() const {
            for (int i = 0; i < N; ++i)
                if (std::isnan(values[i]))
                    return true;
            return false;
        }
        bool hasInf() const {
            for (int i = 0; i < N; ++i)
                if (std::isinf(values[i]))
                    return true;
            return false;
        }
        bool allFinite() const {
            return !hasNaN() && !hasInf();
        }
        bool hasMinus() const {
            for (int i = 0; i < N; ++i)
                if (values[i] < 0)
                    return true;
            return false;
        }
        
        RealType luminance(RGBColorSpace space = RGBColorSpace::sRGB) const {
            RealType sum = 0;
            for (int i = 0; i < N; ++i)
                sum += values[i];
            return sum / N;
        }
        
        // setting "primary" to 1.0 might introduce bias.
        RealType importance(uint16_t selectedLambda) const {
            // I hope a compiler to optimize away this if statement...
            // What I want to do is just only member function specialization of a template class while reusing other function definitions.
            if (N > 1) {
                RealType sum = 0;
                for (int i = 0; i < N; ++i)
                    sum += values[i];
                const RealType primary = 0.9f;
                const RealType marginal = (1 - primary) / (N - 1);
                return sum * marginal + values[selectedLambda] * (primary - marginal);
            }
            else {
                return values[0];
            }
        }
        
        std::string toString() const {
            std::string ret = "(";
            char str[256];
            for (int i = 0; i < N - 1; ++i) {
                sprintf(str, "%g, ", values[i]);
                ret += str;
            }
            sprintf(str, "%g)", values[N - 1]);
            ret += str;
            return ret;
        }
        
        std::string toString(const WavelengthSamplesTemplate<RealType, N> &wls) const {
            std::string ret = "";
            char str[256];
            for (int i = 0; i < N; ++i) {
                sprintf(str, "%g, %g\n", wls[i], values[i]);
                ret += str;
            }
            return ret;
        }

        static const uint32_t NumComponents;
        static const SampledSpectrumTemplate Zero;
        static const SampledSpectrumTemplate One;
        static const SampledSpectrumTemplate Inf;
        static const SampledSpectrumTemplate NaN;
    };
    template <typename RealType, uint32_t N>
    const uint32_t SampledSpectrumTemplate<RealType, N>::NumComponents = N;
    template <typename RealType, uint32_t N>
    const SampledSpectrumTemplate<RealType, N> SampledSpectrumTemplate<RealType, N>::Zero = SampledSpectrumTemplate<RealType, N>(0.0);
    template <typename RealType, uint32_t N>
    const SampledSpectrumTemplate<RealType, N> SampledSpectrumTemplate<RealType, N>::One = SampledSpectrumTemplate<RealType, N>(1.0);
    template <typename RealType, uint32_t N>
    const SampledSpectrumTemplate<RealType, N> SampledSpectrumTemplate<RealType, N>::Inf = SampledSpectrumTemplate<RealType, N>(std::numeric_limits<RealType>::infinity());
    template <typename RealType, uint32_t N>
    const SampledSpectrumTemplate<RealType, N> SampledSpectrumTemplate<RealType, N>::NaN = SampledSpectrumTemplate<RealType, N>(std::numeric_limits<RealType>::quiet_NaN());
    
    template <typename RealType, uint32_t N>
    SLR_API SampledSpectrumTemplate<RealType, N> sqrt(const SampledSpectrumTemplate<RealType, N> &value);
    
    template <typename RealType, uint32_t N>
    SLR_API SampledSpectrumTemplate<RealType, N> exp(const SampledSpectrumTemplate<RealType, N> &value);


    
    template <typename RealType, uint32_t numStrata>
    struct SLR_API DiscretizedSpectrumTemplate {
        RealType values[numStrata];
        
    public:
        DiscretizedSpectrumTemplate(RealType v = 0.0f) { for (int i = 0; i < numStrata; ++i) values[i] = v; }
        DiscretizedSpectrumTemplate(const RealType* vals) { for (int i = 0; i < numStrata; ++i) values[i] = vals[i]; }
        
        DiscretizedSpectrumTemplate operator+() const { return *this; }
        DiscretizedSpectrumTemplate operator-() const {
            RealType vals[numStrata];
            for (int i = 0; i < numStrata; ++i)
                vals[i] = -values[i];
            return DiscretizedSpectrumTemplate(vals);
        }
        
        DiscretizedSpectrumTemplate operator+(const DiscretizedSpectrumTemplate &c) const {
            RealType vals[numStrata];
            for (int i = 0; i < numStrata; ++i)
                vals[i] = values[i] + c.values[i];
            return DiscretizedSpectrumTemplate(vals);
        }
        DiscretizedSpectrumTemplate operator-(const DiscretizedSpectrumTemplate &c) const {
            RealType vals[numStrata];
            for (int i = 0; i < numStrata; ++i)
                vals[i] = values[i] - c.values[i];
            return DiscretizedSpectrumTemplate(vals);
        }
        DiscretizedSpectrumTemplate operator*(const DiscretizedSpectrumTemplate &c) const {
            RealType vals[numStrata];
            for (int i = 0; i < numStrata; ++i)
                vals[i] = values[i] * c.values[i];
            return DiscretizedSpectrumTemplate(vals);
        }
        DiscretizedSpectrumTemplate operator*(RealType s) const {
            RealType vals[numStrata];
            for (int i = 0; i < numStrata; ++i)
                vals[i] = values[i] * s;
            return DiscretizedSpectrumTemplate(vals);
        }
        friend inline DiscretizedSpectrumTemplate operator*(RealType s, const DiscretizedSpectrumTemplate &c) {
            RealType vals[numStrata];
            for (int i = 0; i < numStrata; ++i)
                vals[i] = c.values[i] * s;
            return DiscretizedSpectrumTemplate(vals);
        }
        
        DiscretizedSpectrumTemplate &operator+=(const DiscretizedSpectrumTemplate &c) {
            for (int i = 0; i < numStrata; ++i)
                values[i] += c.values[i];
            return *this;
        }
        DiscretizedSpectrumTemplate &operator*=(const DiscretizedSpectrumTemplate &c) {
            for (int i = 0; i < numStrata; ++i)
                values[i] *= c.values[i];
            return *this;
        }
        DiscretizedSpectrumTemplate &operator*=(RealType s) {
            for (int i = 0; i < numStrata; ++i)
                values[i] *= s;
            return *this;
        }
        
        bool operator==(const DiscretizedSpectrumTemplate &c) const {
            for (int i = 0; i < numStrata; ++i)
                if (values[i] != c.values[i])
                    return false;
            return true;
        }
        bool operator!=(const DiscretizedSpectrumTemplate &c) const {
            for (int i = 0; i < numStrata; ++i)
                if (values[i] != c.values[i])
                    return true;
            return false;
        }
        
        RealType &operator[](unsigned int index) {
            SLRAssert(index < numStrata, "\"index\" is out of range [0, %u].", numStrata - 1);
            return values[index];
        }
        RealType operator[](unsigned int index) const {
            SLRAssert(index < numStrata, "\"index\" is out of range [0, %u].", numStrata - 1);
            return values[index];
        }
        
        RealType maxValue() const {
            RealType maxVal = values[0];
            for (int i = 1; i < numStrata; ++i)
                maxVal = std::fmax(values[i], maxVal);
            return maxVal;
        }
        RealType minValue() const {
            RealType minVal = values[0];
            for (int i = 1; i < numStrata; ++i)
                minVal = std::fmin(values[i], minVal);
            return minVal;
        }
        bool hasNonZero() const {
            for (int i = 0; i < numStrata; ++i)
                if (values[i] != 0)
                    return true;
            return false;
        }
        bool hasNaN() const {
            for (int i = 0; i < numStrata; ++i)
                if (std::isnan(values[i]))
                    return true;
            return false;
        }
        bool hasInf() const {
            for (int i = 0; i < numStrata; ++i)
                if (std::isinf(values[i]))
                    return true;
            return false;
        }
        bool hasMinus() const {
            for (int i = 0; i < numStrata; ++i)
                if (values[i] < 0)
                    return true;
            return false;
        }
        
        RealType luminance(RGBColorSpace space = RGBColorSpace::sRGB) const {
            RealType sum = 0;
            for (int i = 0; i < numStrata; ++i)
                sum += ybar[i] * values[i];
            return sum / integralCMF;
        }
        
        void getRGB(RealType RGB[3], RGBColorSpace space = RGBColorSpace::sRGB) const {
            RealType XYZ[3] = {0, 0, 0};
            for (int i = 0; i < numStrata; ++i) {
                XYZ[0] += xbar[i] * values[i];
                XYZ[1] += ybar[i] * values[i];
                XYZ[2] += zbar[i] * values[i];
            }
            XYZ[0] /= integralCMF;
            XYZ[1] /= integralCMF;
            XYZ[2] /= integralCMF;
            
            switch (space) {
                case RGBColorSpace::sRGB:
                    XYZ_to_sRGB(XYZ, RGB);
                    break;
                default:
                    SLRAssert(false, "Invalid RGB color space is specified.");
                    break;
            }
        }
        
        std::string toString() const {
            std::string ret = "(";
            char str[256];
            for (int i = 0; i < numStrata - 1; ++i) {
                sprintf(str, "%g, ", values[i]);
                ret += str;
            }
            sprintf(str, "%g)", values[numStrata - 1]);
            ret += str;
            return ret;
        }
        
        static const DiscretizedSpectrumTemplate Zero;
        static const DiscretizedSpectrumTemplate One;
        static const DiscretizedSpectrumTemplate Inf;
        static const DiscretizedSpectrumTemplate NaN;
        
        static const uint32_t NumStrata;
        static std::unique_ptr<float[]> xbar;
        static std::unique_ptr<float[]> ybar;
        static std::unique_ptr<float[]> zbar;
        static float integralCMF;
        
        static void init() {
            xbar = std::unique_ptr<float[]>(new float[numStrata]);
            ybar = std::unique_ptr<float[]>(new float[numStrata]);
            zbar = std::unique_ptr<float[]>(new float[numStrata]);
            
            uint32_t sBin = 0;
            float nextP = float(sBin + 1) / numStrata;
            float xSum = 0, xPrev = xbar_2deg[0];
            float ySum = 0, yPrev = ybar_2deg[0];
            float zSum = 0, zPrev = zbar_2deg[0];
            const float interval = 1;// nm
            for (int i = 1; i < NumCMFSamples; ++i) {
                float curP = float(i) / (NumCMFSamples - 1);
                float width = interval;
                float xCur = xbar_2deg[i];
                float yCur = ybar_2deg[i];
                float zCur = zbar_2deg[i];
                if (curP >= nextP) {
                    width = (curP - nextP) * (WavelengthHighBound - WavelengthLowBound);
                    float t = 1 - width / interval;
                    float xIn = xPrev * (1 - t) + xCur * t;
                    float yIn = yPrev * (1 - t) + yCur * t;
                    float zIn = zPrev * (1 - t) + zCur * t;
                    xSum += (xPrev + xIn) * (interval - width) * 0.5f;
                    ySum += (yPrev + yIn) * (interval - width) * 0.5f;
                    zSum += (zPrev + zIn) * (interval - width) * 0.5f;
                    xbar[sBin] = xSum;
                    ybar[sBin] = ySum;
                    zbar[sBin] = zSum;
                    
                    xSum = ySum = zSum = 0;
                    xPrev = xIn;
                    yPrev = yIn;
                    zPrev = zIn;
                    
                    ++sBin;
                    nextP = float(sBin + 1) / numStrata;
                }
                xSum += (xPrev + xCur) * width * 0.5f;
                ySum += (yPrev + yCur) * width * 0.5f;
                zSum += (zPrev + zCur) * width * 0.5f;
                xPrev = xCur;
                yPrev = yCur;
                zPrev = zCur;
            }
            
            integralCMF = 0.0f;
            for (int i = 0; i < numStrata; ++i)
                integralCMF += ybar[i];
        }
    };
    template <typename RealType, uint32_t numStrata>
    const uint32_t DiscretizedSpectrumTemplate<RealType, numStrata>::NumStrata = numStrata;
    template <typename RealType, uint32_t numStrata>
    std::unique_ptr<float[]> DiscretizedSpectrumTemplate<RealType, numStrata>::xbar;
    template <typename RealType, uint32_t numStrata>
    std::unique_ptr<float[]> DiscretizedSpectrumTemplate<RealType, numStrata>::ybar;
    template <typename RealType, uint32_t numStrata>
    std::unique_ptr<float[]> DiscretizedSpectrumTemplate<RealType, numStrata>::zbar;
    template <typename RealType, uint32_t numStrata>
    float DiscretizedSpectrumTemplate<RealType, numStrata>::integralCMF;
    template <typename RealType, uint32_t numStrata>
    const DiscretizedSpectrumTemplate<RealType, numStrata> DiscretizedSpectrumTemplate<RealType, numStrata>::Zero = DiscretizedSpectrumTemplate<RealType, numStrata>(0.0);
    template <typename RealType, uint32_t numStrata>
    const DiscretizedSpectrumTemplate<RealType, numStrata> DiscretizedSpectrumTemplate<RealType, numStrata>::One = DiscretizedSpectrumTemplate<RealType, numStrata>(1.0);
    template <typename RealType, uint32_t numStrata>
    const DiscretizedSpectrumTemplate<RealType, numStrata> DiscretizedSpectrumTemplate<RealType, numStrata>::Inf = DiscretizedSpectrumTemplate<RealType, numStrata>(std::numeric_limits<RealType>::infinity());
    template <typename RealType, uint32_t numStrata>
    const DiscretizedSpectrumTemplate<RealType, numStrata> DiscretizedSpectrumTemplate<RealType, numStrata>::NaN = DiscretizedSpectrumTemplate<RealType, numStrata>(std::numeric_limits<RealType>::quiet_NaN());


    
    template <typename RealType, uint32_t numStrata>
    class SLR_API SpectrumStorageTemplate {
        typedef DiscretizedSpectrumTemplate<RealType, numStrata> ValueType;
        CompensatedSum<ValueType> value;
        
    public:
        SpectrumStorageTemplate(const ValueType &v = ValueType::Zero) :
        value(v) {}
        
        template <uint32_t N>
        SpectrumStorageTemplate &add(const WavelengthSamplesTemplate<RealType, N> &wls, const SampledSpectrumTemplate<RealType, N> &val) {
            const RealType recBinWidth = numStrata / (WavelengthHighBound - WavelengthLowBound);
            ValueType addend(0.0);
            for (int i = 0; i < WavelengthSamplesTemplate<RealType, N>::NumComponents; ++i) {
                uint32_t sBin = std::min(uint32_t((wls[i] - WavelengthLowBound) / (WavelengthHighBound - WavelengthLowBound) * numStrata), numStrata - 1);
                addend[sBin] += val[i] * recBinWidth;
            }
            value += addend;
            return *this;
        }
        
        CompensatedSum<ValueType> &getValue() {
            return value;
        }
    };
}

#endif
