//
//  spectrum_types.h
//
//  Created by 渡部 心 on 2015/09/13.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_spectrum_types__
#define __SLR_spectrum_types__

#include "../defines.h"
#include "../declarations.h"
#include "spectrum_base.h"
#include "rgb_types.h"
#include "CompensatedSum.h"

namespace SLR {
    template <typename RealType, uint32_t NumSpectralSamples>
    struct SLR_API WavelengthSamplesTemplate {
        enum Flag : uint16_t {
            WavelengthIsSelected = 0x01,
        };

        RealType lambdas[NumSpectralSamples];
        uint16_t selectedLambdaIndex;
        uint16_t flags;
        
        WavelengthSamplesTemplate() : selectedLambdaIndex(0), flags(0) {};
        WavelengthSamplesTemplate(const RealType* values) {
            for (int i = 0; i < NumSpectralSamples; ++i)
                lambdas[i] = values[i];
            selectedLambdaIndex = 0;
            flags = 0;
        }
        WavelengthSamplesTemplate(const WavelengthSamplesTemplate &wls) {
            for (int i = 0; i < NumSpectralSamples; ++i)
                lambdas[i] = wls.lambdas[i];
            selectedLambdaIndex = wls.selectedLambdaIndex;
            flags = wls.flags;
        }
        
        RealType &operator[](uint32_t index) {
            SLRAssert(index < NumSpectralSamples, "\"index\" is out of range [0, %u].", NumSpectralSamples - 1);
            return lambdas[index];
        }
        RealType operator[](uint32_t index) const {
            SLRAssert(index < NumSpectralSamples, "\"index\" is out of range [0, %u].", NumSpectralSamples - 1);
            return lambdas[index];
        }
        
        bool wavelengthSelected() const {
            return (flags & WavelengthIsSelected) != 0;
        }
        
        RealType selectedWavelength() const {
            return lambdas[selectedLambdaIndex];
        }
        
        static WavelengthSamplesTemplate createWithEqualOffsets(RealType offset, RealType uLambda, RealType* PDF) {
            SLRAssert(offset >= 0 && offset < 1, "\"offset\" must be in range [0, 1).");
            SLRAssert(uLambda >= 0 && uLambda < 1, "\"uLambda\" must be in range [0, 1).");
            WavelengthSamplesTemplate wls;
            for (int i = 0; i < NumSpectralSamples; ++i)
                wls.lambdas[i] = WavelengthLowBound + (WavelengthHighBound - WavelengthLowBound) * (i + offset) / NumSpectralSamples;
            wls.selectedLambdaIndex = std::min(uint16_t(NumSpectralSamples * uLambda), uint16_t(NumSpectralSamples - 1));
            wls.flags = 0;
            *PDF = NumSpectralSamples / (WavelengthHighBound - WavelengthLowBound);
            return wls;
        }
        
        static const uint32_t NumComponents;
    };

    
    
    template <typename RealType, uint32_t NumSpectralSamples>
    class SLR_API ContinuousSpectrumTemplate {
    public:
        virtual ~ContinuousSpectrumTemplate() { }
        
        virtual void calcBounds(uint32_t numBins, RealType* bounds) const = 0;
        virtual SampledSpectrumTemplate<RealType, NumSpectralSamples> evaluate(const WavelengthSamplesTemplate<RealType, NumSpectralSamples> &wls) const = 0;
        virtual void evaluate(const RealType* wavelengths, uint32_t numSamples, RealType* values) const = 0;
        virtual void convertToXYZ(RealType XYZ[3]) const = 0;
        virtual ContinuousSpectrumTemplate* createScaledAndOffset(RealType scale, RealType offset) const = 0;
        
        RGBTemplate<RealType> convertToRGB(SpectrumType spType) const;
    };
    
    
    
    template <typename RealType, uint32_t NumSpectralSamples>
    class SLR_API RegularContinuousSpectrumTemplate : public ContinuousSpectrumTemplate<RealType, NumSpectralSamples> {
        RealType m_minLambda, m_maxLambda;
        uint32_t m_numSamples;
        RealType* m_values;
        
    public:
        RegularContinuousSpectrumTemplate(RealType minWL, RealType maxWL, const RealType* vals, uint32_t numVals) : m_minLambda(minWL), m_maxLambda(maxWL), m_numSamples(numVals) {
            m_values = new RealType[m_numSamples];
            for (int i = 0; i < m_numSamples; ++i)
                m_values[i] = vals[i];
        }
        ~RegularContinuousSpectrumTemplate() {
            delete[] m_values;
        }
        
        void calcBounds(uint32_t numBins, RealType* bounds) const override;
        SampledSpectrumTemplate<RealType, NumSpectralSamples> evaluate(const WavelengthSamplesTemplate<RealType, NumSpectralSamples> &wls) const override;
        void evaluate(const RealType* wavelengths, uint32_t numSamples, RealType* values) const override;
        void convertToXYZ(RealType XYZ[3]) const override;
        ContinuousSpectrumTemplate<RealType, NumSpectralSamples>* createScaledAndOffset(RealType scale, RealType offset) const override;
    };

    
    
    template <typename RealType, uint32_t NumSpectralSamples>
    class SLR_API IrregularContinuousSpectrumTemplate : public ContinuousSpectrumTemplate<RealType, NumSpectralSamples> {
        uint32_t m_numSamples;
        RealType* m_lambdas;
        RealType* m_values;
        
    public:
        IrregularContinuousSpectrumTemplate(const RealType* wls, const RealType* vals, uint32_t numVals) : m_numSamples(numVals) {
            m_lambdas = new RealType[m_numSamples];
            m_values = new RealType[m_numSamples];
            for (int i = 0; i < m_numSamples; ++i) {
                m_lambdas[i] = wls[i];
                m_values[i] = vals[i];
            }
        }
        ~IrregularContinuousSpectrumTemplate() {
            delete[] m_lambdas;
            delete[] m_values;
        }
        
        void calcBounds(uint32_t numBins, RealType* bounds) const override;
        SampledSpectrumTemplate<RealType, NumSpectralSamples> evaluate(const WavelengthSamplesTemplate<RealType, NumSpectralSamples> &wls) const override;
        void evaluate(const RealType* wavelengths, uint32_t numSamples, RealType* values) const override;
        void convertToXYZ(RealType XYZ[3]) const override;
        ContinuousSpectrumTemplate<RealType, NumSpectralSamples>* createScaledAndOffset(RealType scale, RealType offset) const override;
    };

    
    
    // JP: これは理想的には不要である。
    //     環境テクスチャーをuvs16Fx3(uvsA16Fx4)フォーマットに変換して保持する際に、典型的な値の場合、変換後のsの値が容易にhalf floatの限界に達するため、
    //     適当なスケール値をかけて小さな値にする。
    // EN: This is not ideally needed.
    //     When converting an environment texture into uvs16Fx3 (uvsA16Fx4) format and holding it, 
    //     a resulting s value from a typical value easily reaches the limit of half float therefore make it small by multiplying an appropriate scaling value.  
#define UPSAMPLED_CONTINOUS_SPECTRUM_SCALE_FACTOR (0.009355121400914532) // corresponds to cancel dividing by EqualEnergyReflectance.
    
    // References
    // Physically Meaningful Rendering using Tristimulus Colours
    template <typename RealType, uint32_t NumSpectralSamples>
    class SLR_API UpsampledContinuousSpectrumTemplate : public ContinuousSpectrumTemplate<RealType, NumSpectralSamples> {
        // Grid cells. Laid out in row-major format.
        // num_points = 0 for cells without data points.
        struct spectrum_grid_cell_t {
            uint8_t inside;
            uint8_t num_points;
            uint8_t idx[6];
        };
        
        // Grid data points.
        struct spectrum_data_point_t {
            RealType xystar[2];
            RealType uv[2];
            RealType spectrum[95]; // X+Y+Z = 1
        };
        
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
        
        void calcBounds(uint32_t numBins, RealType* bounds) const override;
        SampledSpectrumTemplate<RealType, NumSpectralSamples> evaluate(const WavelengthSamplesTemplate<RealType, NumSpectralSamples> &wls) const override;
        void evaluate(const RealType* wavelengths, uint32_t numSamples, RealType* values) const override;
        void convertToXYZ(RealType XYZ[3]) const override;
        ContinuousSpectrumTemplate<RealType, NumSpectralSamples>* createScaledAndOffset(RealType scale, RealType offset) const override;
        
        static const RealType MinWavelength;
        static const RealType MaxWavelength;
        static const uint32_t NumWavelengthSamples;
        static const uint32_t GridWidth;
        static const uint32_t GridHeight;
        
        // This is 1 over the integral over either CMF.
        // Spectra can be mapped so that xyz=(1,1,1) is converted to constant 1 by
        // dividing by this value. This is important for valid reflectances.
        static const RealType EqualEnergyReflectance;
        static const spectrum_grid_cell_t spectrum_grid[];
        static const spectrum_data_point_t spectrum_data_points[];
        
        static inline void xy_to_uv(const RealType xy[2], RealType uv[2]) {
            uv[0] = 16.730260708356887 * xy[0] + 7.7801960340706 * xy[1] - 2.170152247475828,
            uv[1] = -7.530081094743006 * xy[0] + 16.192422314095225 * xy[1] + 1.1125529268825947;
        }
        
        static inline void uv_to_xy(const RealType uv[2], RealType xy[2]) {
            xy[0] = 0.0491440520940413 * uv[0] - 0.02361291916573777 * uv[1] + 0.13292069743203658;
            xy[1] = 0.022853819546830627 * uv[0] + 0.05077639329371236 * uv[1] - 0.006895157122499944;
        }
        
        static inline void sRGB_to_uvs(SpectrumType spType, const RealType rgb[3], RealType uvs[3]) {
            RealType xyz[3];
            switch (spType) {
                case SpectrumType::Reflectance:
                    sRGB_E_to_XYZ(rgb, xyz);
                    break;
                case SpectrumType::LightSource:
                    sRGB_to_XYZ(rgb, xyz);
                    break;
                case SpectrumType::IndexOfRefraction:
                    sRGB_E_to_XYZ(rgb, xyz);
                    break;
                default:
                    break;
            }
            RealType xy[2];
            RealType b = xyz[0] + xyz[1] + xyz[2];
            xy[0] = xyz[0] / b;
            xy[1] = xyz[1] / b;
            if (b == 0)
                xy[0] = xy[1] = 1.0f / 3.0;
            xy_to_uv(xy, uvs);
            uvs[2] = b / EqualEnergyReflectance;
        }
        
        static inline void uvs_to_sRGB(SpectrumType spType, const RealType uvs[3], RealType rgb[3]) {
            RealType xy[2];
            uv_to_xy(uvs, xy);
            RealType b = uvs[2] * EqualEnergyReflectance;
            RealType XYZ[3];
            XYZ[0] = xy[0] * b;
            XYZ[1] = xy[1] * b;
            XYZ[2] = b - XYZ[0] - XYZ[1];
            switch (spType) {
                case SpectrumType::Reflectance:
                    XYZ_to_sRGB_E(XYZ, rgb);
                    break;
                case SpectrumType::LightSource:
                    XYZ_to_sRGB(XYZ, rgb);
                    break;
                case SpectrumType::IndexOfRefraction:
                    XYZ_to_sRGB_E(XYZ, rgb);
                    break;
                default:
                    break;
            }
        }
        
        static inline RealType uvs_to_luminance(const RealType uvs[3]) {
            RealType xy[2];
            uv_to_xy(uvs, xy);
            RealType b = uvs[2] * EqualEnergyReflectance;
            return xy[1] * b;
        }
    };

    template <typename RealType, uint32_t NumSpectralSamples>
    const RealType UpsampledContinuousSpectrumTemplate<RealType, NumSpectralSamples>::EqualEnergyReflectance = 0.009355121400914532;
    
    
    
    template <typename RealType, uint32_t NumSpectralSamples>
    class SLR_API ScaledAndOffsetUpsampledContinuousSpectrumTemplate : public ContinuousSpectrumTemplate<RealType, NumSpectralSamples> {
        const UpsampledContinuousSpectrumTemplate<RealType, NumSpectralSamples> m_baseSpectrum;
        RealType m_scale;
        RealType m_offset;
        
    public:
        ScaledAndOffsetUpsampledContinuousSpectrumTemplate(const UpsampledContinuousSpectrumTemplate<RealType, NumSpectralSamples> &baseSpectrum, 
                                                           RealType scale, RealType offset) :
        m_baseSpectrum(baseSpectrum), m_scale(scale), m_offset(offset) { }
        
        void calcBounds(uint32_t numBins, RealType* bounds) const override;
        SampledSpectrumTemplate<RealType, NumSpectralSamples> evaluate(const WavelengthSamplesTemplate<RealType, NumSpectralSamples> &wls) const override;
        void evaluate(const RealType* wavelengths, uint32_t numSamples, RealType* values) const override;
        void convertToXYZ(RealType XYZ[3]) const override;
        ContinuousSpectrumTemplate<RealType, NumSpectralSamples>* createScaledAndOffset(RealType scale, RealType offset) const override;
    };

    

    template <typename RealType, uint32_t NumSpectralSamples>
    struct SLR_API SampledSpectrumTemplate {
        RealType values[NumSpectralSamples];

        SampledSpectrumTemplate(RealType v = 0.0f) { for (int i = 0; i < NumSpectralSamples; ++i) values[i] = v; }
        SampledSpectrumTemplate(const RealType* vals) { for (int i = 0; i < NumSpectralSamples; ++i) values[i] = vals[i]; }
        
        SampledSpectrumTemplate operator+() const { return *this; };
        SampledSpectrumTemplate operator-() const {
            RealType vals[NumSpectralSamples];
            for (int i = 0; i < NumSpectralSamples; ++i)
                vals[i] = -values[i];
            return SampledSpectrumTemplate(vals);
        }
        
        SampledSpectrumTemplate operator+(const SampledSpectrumTemplate &c) const {
            RealType vals[NumSpectralSamples];
            for (int i = 0; i < NumSpectralSamples; ++i)
                vals[i] = values[i] + c.values[i];
            return SampledSpectrumTemplate(vals);
        }
        SampledSpectrumTemplate operator-(const SampledSpectrumTemplate &c) const {
            RealType vals[NumSpectralSamples];
            for (int i = 0; i < NumSpectralSamples; ++i)
                vals[i] = values[i] - c.values[i];
            return SampledSpectrumTemplate(vals);
        }
        SampledSpectrumTemplate operator*(const SampledSpectrumTemplate &c) const {
            RealType vals[NumSpectralSamples];
            for (int i = 0; i < NumSpectralSamples; ++i)
                vals[i] = values[i] * c.values[i];
            return SampledSpectrumTemplate(vals);
        }
        SampledSpectrumTemplate operator/(const SampledSpectrumTemplate &c) const {
            RealType vals[NumSpectralSamples];
            for (int i = 0; i < NumSpectralSamples; ++i)
                vals[i] = values[i] / c.values[i];
            return SampledSpectrumTemplate(vals);
        }
        SampledSpectrumTemplate safeDivide(const SampledSpectrumTemplate &c) const {
            RealType vals[NumSpectralSamples];
            for (int i = 0; i < NumSpectralSamples; ++i)
                vals[i] = c.values[i] > 0 ? values[i] / c.values[i] : 0.0f;
            return SampledSpectrumTemplate(vals);
        }
        SampledSpectrumTemplate operator*(RealType s) const {
            RealType vals[NumSpectralSamples];
            for (int i = 0; i < NumSpectralSamples; ++i)
                vals[i] = values[i] * s;
            return SampledSpectrumTemplate(vals);
        }
        SampledSpectrumTemplate operator/(RealType s) const {
            RealType vals[NumSpectralSamples];
            RealType r = 1 / s;
            for (int i = 0; i < NumSpectralSamples; ++i)
                vals[i] = values[i] * r;
            return SampledSpectrumTemplate(vals);
        }
        friend inline SampledSpectrumTemplate operator*(RealType s, const SampledSpectrumTemplate &c) {
            RealType vals[NumSpectralSamples];
            for (int i = 0; i < NumSpectralSamples; ++i)
                vals[i] = c.values[i] * s;
            return SampledSpectrumTemplate(vals);
        }
        
        SampledSpectrumTemplate &operator+=(const SampledSpectrumTemplate &c) {
            for (int i = 0; i < NumSpectralSamples; ++i)
                values[i] += c.values[i];
            return *this;
        }
        SampledSpectrumTemplate &operator-=(const SampledSpectrumTemplate &c) {
            for (int i = 0; i < NumSpectralSamples; ++i)
                values[i] -= c.values[i];
            return *this;
        }
        SampledSpectrumTemplate &operator*=(const SampledSpectrumTemplate &c) {
            for (int i = 0; i < NumSpectralSamples; ++i)
                values[i] *= c.values[i];
            return *this;
        }
        SampledSpectrumTemplate &operator/=(const SampledSpectrumTemplate &c) {
            for (int i = 0; i < NumSpectralSamples; ++i)
                values[i] /= c.values[i];
            return *this;
        }
        SampledSpectrumTemplate &operator*=(RealType s) {
            for (int i = 0; i < NumSpectralSamples; ++i)
                values[i] *= s;
            return *this;
        }
        SampledSpectrumTemplate &operator/=(RealType s) {
            RealType r = 1 / s;
            for (int i = 0; i < NumSpectralSamples; ++i)
                values[i] *= r;
            return *this;
        }
        
        bool operator==(const SampledSpectrumTemplate &c) const {
            for (int i = 0; i < NumSpectralSamples; ++i)
                if (values[i] != c.values[i])
                    return false;
            return true;
        }
        bool operator!=(const SampledSpectrumTemplate &c) const {
            for (int i = 0; i < NumSpectralSamples; ++i)
                if (values[i] != c.values[i])
                    return true;
            return false;
        }
        
        RealType &operator[](unsigned int index) {
            SLRAssert(index < NumSpectralSamples, "\"index\" is out of range [0, %u].", NumSpectralSamples - 1);
            return values[index];
        }
        RealType operator[](unsigned int index) const {
            SLRAssert(index < NumSpectralSamples, "\"index\" is out of range [0, %u].", NumSpectralSamples - 1);
            return values[index];
        }
        
        RealType avgValue() const {
            RealType sumVal = values[0];
            for (int i = 1; i < NumSpectralSamples; ++i)
                sumVal += values[i];
            return sumVal / NumSpectralSamples;
        }
        RealType maxValue() const {
            RealType maxVal = values[0];
            for (int i = 1; i < NumSpectralSamples; ++i)
                maxVal = std::fmax(values[i], maxVal);
            return maxVal;
        }
        RealType minValue() const {
            RealType minVal = values[0];
            for (int i = 1; i < NumSpectralSamples; ++i)
                minVal = std::fmin(values[i], minVal);
            return minVal;
        }
        bool hasNonZero() const {
            for (int i = 0; i < NumSpectralSamples; ++i)
                if (values[i] != 0)
                    return true;
            return false;
        }
        bool hasNaN() const {
            for (int i = 0; i < NumSpectralSamples; ++i)
                if (std::isnan(values[i]))
                    return true;
            return false;
        }
        bool hasInf() const {
            for (int i = 0; i < NumSpectralSamples; ++i)
                if (std::isinf(values[i]))
                    return true;
            return false;
        }
        bool allFinite() const {
            return !hasNaN() && !hasInf();
        }
        bool hasNegative() const {
            for (int i = 0; i < NumSpectralSamples; ++i)
                if (values[i] < 0)
                    return true;
            return false;
        }
        
        RealType luminance(RGBColorSpace space = RGBColorSpace::sRGB) const {
            RealType sum = 0;
            for (int i = 0; i < NumSpectralSamples; ++i)
                sum += values[i];
            return sum / NumSpectralSamples;
        }
        
        // setting "primary" to 1.0 might introduce bias.
        RealType importance(uint16_t selectedLambda) const {
            // I hope a compiler to optimize away this if statement...
            // What I want to do is just only member function specialization of a template class while reusing other function definitions.
            if (NumSpectralSamples > 1) {
                RealType sum = 0;
                for (int i = 0; i < NumSpectralSamples; ++i)
                    sum += values[i];
                const RealType primary = 0.9f;
                const RealType marginal = (1 - primary) / (NumSpectralSamples - 1);
                return sum * marginal + values[selectedLambda] * (primary - marginal);
            }
            else {
                return values[0];
            }
        }
        
        std::string toString() const {
            std::string ret = "(";
            char str[256];
            for (int i = 0; i < NumSpectralSamples - 1; ++i) {
                sprintf(str, "%g, ", values[i]);
                ret += str;
            }
            sprintf(str, "%g)", values[NumSpectralSamples - 1]);
            ret += str;
            return ret;
        }
        
        std::string toString(const WavelengthSamplesTemplate<RealType, NumSpectralSamples> &wls) const {
            std::string ret = "";
            char str[256];
            for (int i = 0; i < NumSpectralSamples; ++i) {
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
    
    template <typename RealType, uint32_t NumSpectralSamples>
    SLR_API SampledSpectrumTemplate<RealType, NumSpectralSamples> min(const SampledSpectrumTemplate<RealType, NumSpectralSamples> &value, RealType minValue);
    
    template <typename RealType, uint32_t NumSpectralSamples>
    SLR_API SampledSpectrumTemplate<RealType, NumSpectralSamples> max(const SampledSpectrumTemplate<RealType, NumSpectralSamples> &value, RealType maxValue);
    
    template <typename RealType, uint32_t NumSpectralSamples>
    SLR_API SampledSpectrumTemplate<RealType, NumSpectralSamples> lerp(const SampledSpectrumTemplate<RealType, NumSpectralSamples> &v0, 
                                                                       const SampledSpectrumTemplate<RealType, NumSpectralSamples> &v1, 
                                                                       RealType t);
    
    template <typename RealType, uint32_t NumSpectralSamples>
    SLR_API SampledSpectrumTemplate<RealType, NumSpectralSamples> sqrt(const SampledSpectrumTemplate<RealType, NumSpectralSamples> &value);
    
    template <typename RealType, uint32_t NumSpectralSamples>
    SLR_API SampledSpectrumTemplate<RealType, NumSpectralSamples> exp(const SampledSpectrumTemplate<RealType, NumSpectralSamples> &value);


    
    template <typename RealType, uint32_t NumStrataForStorage>
    struct SLR_API DiscretizedSpectrumTemplate {
        RealType values[NumStrataForStorage];
        
    public:
        DiscretizedSpectrumTemplate(RealType v = 0.0f) { for (int i = 0; i < NumStrataForStorage; ++i) values[i] = v; }
        DiscretizedSpectrumTemplate(const RealType* vals) { for (int i = 0; i < NumStrataForStorage; ++i) values[i] = vals[i]; }
        
        DiscretizedSpectrumTemplate operator+() const { return *this; }
        DiscretizedSpectrumTemplate operator-() const {
            RealType vals[NumStrataForStorage];
            for (int i = 0; i < NumStrataForStorage; ++i)
                vals[i] = -values[i];
            return DiscretizedSpectrumTemplate(vals);
        }
        
        DiscretizedSpectrumTemplate operator+(const DiscretizedSpectrumTemplate &c) const {
            RealType vals[NumStrataForStorage];
            for (int i = 0; i < NumStrataForStorage; ++i)
                vals[i] = values[i] + c.values[i];
            return DiscretizedSpectrumTemplate(vals);
        }
        DiscretizedSpectrumTemplate operator-(const DiscretizedSpectrumTemplate &c) const {
            RealType vals[NumStrataForStorage];
            for (int i = 0; i < NumStrataForStorage; ++i)
                vals[i] = values[i] - c.values[i];
            return DiscretizedSpectrumTemplate(vals);
        }
        DiscretizedSpectrumTemplate operator*(const DiscretizedSpectrumTemplate &c) const {
            RealType vals[NumStrataForStorage];
            for (int i = 0; i < NumStrataForStorage; ++i)
                vals[i] = values[i] * c.values[i];
            return DiscretizedSpectrumTemplate(vals);
        }
        DiscretizedSpectrumTemplate operator*(RealType s) const {
            RealType vals[NumStrataForStorage];
            for (int i = 0; i < NumStrataForStorage; ++i)
                vals[i] = values[i] * s;
            return DiscretizedSpectrumTemplate(vals);
        }
        friend inline DiscretizedSpectrumTemplate operator*(RealType s, const DiscretizedSpectrumTemplate &c) {
            RealType vals[NumStrataForStorage];
            for (int i = 0; i < NumStrataForStorage; ++i)
                vals[i] = c.values[i] * s;
            return DiscretizedSpectrumTemplate(vals);
        }
        
        DiscretizedSpectrumTemplate &operator+=(const DiscretizedSpectrumTemplate &c) {
            for (int i = 0; i < NumStrataForStorage; ++i)
                values[i] += c.values[i];
            return *this;
        }
        DiscretizedSpectrumTemplate &operator*=(const DiscretizedSpectrumTemplate &c) {
            for (int i = 0; i < NumStrataForStorage; ++i)
                values[i] *= c.values[i];
            return *this;
        }
        DiscretizedSpectrumTemplate &operator*=(RealType s) {
            for (int i = 0; i < NumStrataForStorage; ++i)
                values[i] *= s;
            return *this;
        }
        
        bool operator==(const DiscretizedSpectrumTemplate &c) const {
            for (int i = 0; i < NumStrataForStorage; ++i)
                if (values[i] != c.values[i])
                    return false;
            return true;
        }
        bool operator!=(const DiscretizedSpectrumTemplate &c) const {
            for (int i = 0; i < NumStrataForStorage; ++i)
                if (values[i] != c.values[i])
                    return true;
            return false;
        }
        
        RealType &operator[](unsigned int index) {
            SLRAssert(index < NumStrataForStorage, "\"index\" is out of range [0, %u].", NumStrataForStorage - 1);
            return values[index];
        }
        RealType operator[](unsigned int index) const {
            SLRAssert(index < NumStrataForStorage, "\"index\" is out of range [0, %u].", NumStrataForStorage - 1);
            return values[index];
        }
        
        RealType maxValue() const {
            RealType maxVal = values[0];
            for (int i = 1; i < NumStrataForStorage; ++i)
                maxVal = std::fmax(values[i], maxVal);
            return maxVal;
        }
        RealType minValue() const {
            RealType minVal = values[0];
            for (int i = 1; i < NumStrataForStorage; ++i)
                minVal = std::fmin(values[i], minVal);
            return minVal;
        }
        bool hasNonZero() const {
            for (int i = 0; i < NumStrataForStorage; ++i)
                if (values[i] != 0)
                    return true;
            return false;
        }
        bool hasNaN() const {
            for (int i = 0; i < NumStrataForStorage; ++i)
                if (std::isnan(values[i]))
                    return true;
            return false;
        }
        bool hasInf() const {
            for (int i = 0; i < NumStrataForStorage; ++i)
                if (std::isinf(values[i]))
                    return true;
            return false;
        }
        bool hasNegative() const {
            for (int i = 0; i < NumStrataForStorage; ++i)
                if (values[i] < 0)
                    return true;
            return false;
        }
        
        RealType luminance(RGBColorSpace space = RGBColorSpace::sRGB) const {
            RealType sum = 0;
            for (int i = 0; i < NumStrataForStorage; ++i)
                sum += ybar[i] * values[i];
            return sum / integralCMF;
        }
        
        void getRGB(RealType RGB[3], RGBColorSpace space = RGBColorSpace::sRGB) const {
            RealType XYZ[3] = {0, 0, 0};
            for (int i = 0; i < NumStrataForStorage; ++i) {
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
            for (int i = 0; i < NumStrataForStorage - 1; ++i) {
                sprintf(str, "%g, ", values[i]);
                ret += str;
            }
            sprintf(str, "%g)", values[NumStrataForStorage - 1]);
            ret += str;
            return ret;
        }
        
        static const DiscretizedSpectrumTemplate Zero;
        static const DiscretizedSpectrumTemplate One;
        static const DiscretizedSpectrumTemplate Inf;
        static const DiscretizedSpectrumTemplate NaN;
        
        static const uint32_t NumStrata;
        static std::unique_ptr<RealType[]> xbar;
        static std::unique_ptr<RealType[]> ybar;
        static std::unique_ptr<RealType[]> zbar;
        static RealType integralCMF;
        
        static void initialize();
    };


    
    template <typename RealType, uint32_t NumStrataForStorage>
    class SLR_API SpectrumStorageTemplate {
        typedef DiscretizedSpectrumTemplate<RealType, NumStrataForStorage> ValueType;
        CompensatedSum<ValueType> value;
        
    public:
        SpectrumStorageTemplate(const ValueType &v = ValueType::Zero) :
        value(v) {}
        
        template <uint32_t N>
        SpectrumStorageTemplate &add(const WavelengthSamplesTemplate<RealType, N> &wls, const SampledSpectrumTemplate<RealType, N> &val) {
            const RealType recBinWidth = NumStrataForStorage / (WavelengthHighBound - WavelengthLowBound);
            ValueType addend(0.0);
            for (int i = 0; i < WavelengthSamplesTemplate<RealType, N>::NumComponents; ++i) {
                uint32_t sBin = std::min(uint32_t((wls[i] - WavelengthLowBound) / (WavelengthHighBound - WavelengthLowBound) * NumStrataForStorage), NumStrataForStorage - 1);
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

#endif /* __SLR_spectrum_types__ */
