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
            LambdaSelected = 0x01,
        };
        RealType lambdas[N];
        uint16_t selectedLambda;
        uint16_t flags;
        static const uint32_t NumComponents;
        
        WavelengthSamplesTemplate() : selectedLambda(0), flags(0) {};
        WavelengthSamplesTemplate(const RealType* values) {
            for (int i = 0; i < N; ++i)
                lambdas[i] = values[i];
            flags = 0;
        };
        
        RealType &operator[](uint32_t index) {
            SLRAssert(index < N, "\"index\" is out of range [0, %u].", N - 1);
            return lambdas[index];
        };
        RealType operator[](uint32_t index) const {
            SLRAssert(index < N, "\"index\" is out of range [0, %u].", N - 1);
            return lambdas[index];
        };
        
        bool lambdaSelected() const {
            return (flags & LambdaSelected) != 0;
        };
        
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
        };
    };
    template <typename RealType, uint32_t N>
    const uint32_t WavelengthSamplesTemplate<RealType, N>::NumComponents = N;

    
    template <typename RealType, uint32_t N>
    struct SLR_API ContinuousSpectrumTemplate {
        virtual SampledSpectrumTemplate<RealType, N> evaluate(const WavelengthSamplesTemplate<RealType, N> &wls) const = 0;
    };
    
    template <typename RealType, uint32_t N>
    struct SLR_API RegularContinuousSpectrumTemplate : public ContinuousSpectrumTemplate<RealType, N> {
        RealType minLambda, maxLambda;
        uint32_t numSamples;
        RealType* values;
        
        RegularContinuousSpectrumTemplate(RealType minWL, RealType maxWL, const RealType* vals, uint32_t numVals) : minLambda(minWL), maxLambda(maxWL), numSamples(numVals) {
            values = new RealType[numSamples];
            for (int i = 0; i < numSamples; ++i)
                values[i] = vals[i];
        };
        ~RegularContinuousSpectrumTemplate() {
            delete[] values;
        };
        
        SampledSpectrumTemplate<RealType, N> evaluate(const WavelengthSamplesTemplate<RealType, N> &wls) const override {
            SampledSpectrumTemplate<RealType, N> ret(0.0f);
            for (int i = 0; i < WavelengthSamplesTemplate<RealType, N>::NumComponents; ++i) {
                RealType binF = (wls[i] - minLambda) / (maxLambda - minLambda) * (numSamples - 1);
                if (binF <= 0.0f) {
                    ret[i] = values[0];
                    continue;
                }
                else if (binF >= numSamples - 1) {
                    ret[i] = values[numSamples - 1];
                    continue;
                }
                int32_t bin = int32_t(binF);
                SLRAssert(bin >= 0 && bin < numSamples - 1, "invalid bin index.");
                RealType t = binF - bin;
                ret[i] = (1 - t) * values[bin] + t * values[bin + 1];
            }
            return ret;
        };
    };

    template <typename RealType, uint32_t N>
    struct SLR_API IrregularContinuousSpectrumTemplate : public ContinuousSpectrumTemplate<RealType, N> {
        uint32_t numSamples;
        RealType* lambdas;
        RealType* values;
        
        IrregularContinuousSpectrumTemplate(const RealType* wls, const RealType* vals, uint32_t numVals) : numSamples(numVals) {
            lambdas = new RealType[numSamples];
            values = new RealType[numSamples];
            for (int i = 0; i < numSamples; ++i) {
                lambdas[i] = wls[i];
                values[i] = vals[i];
            }
        };
        ~IrregularContinuousSpectrumTemplate() {
            delete[] lambdas;
            delete[] values;
        };
        
        SampledSpectrumTemplate<RealType, N> evaluate(const WavelengthSamplesTemplate<RealType, N> &wls) const override {
            SampledSpectrumTemplate<RealType, N> ret(0.0f);
            uint32_t searchBase = 0;
            for (int i = 0; i < WavelengthSamplesTemplate<RealType, N>::NumComponents; ++i) {
                int32_t lowIdx = std::max((int32_t)std::distance(lambdas, std::lower_bound(lambdas + searchBase, lambdas + numSamples, wls[i])) - 1, 0);
                searchBase = lowIdx;
                if (lowIdx >= numSamples - 1) {
                    ret[i] = values[numSamples - 1];
                    continue;
                }
                RealType t = (wls[i] - lambdas[lowIdx]) / (lambdas[lowIdx + 1] - lambdas[lowIdx]);
                if (t <= 0.0f) {
                    ret[i] = values[0];
                    continue;
                }
                SLRAssert(t >= 0 && t <= 1, "invalid interpolation coefficient.");
                ret[i] = (1 - t) * values[lowIdx] + t * values[lowIdx + 1];
            }
            return ret;
        };
    };

    template <typename RealType, uint32_t N>
    struct SLR_API UpsampledContinuousSpectrumTemplate : public ContinuousSpectrumTemplate<RealType, N> {
        RealType u, v, scale;
        
        UpsampledContinuousSpectrumTemplate(RealType uu, RealType vv, RealType ss) : u(uu), v(vv), scale(ss) {};
        
        UpsampledContinuousSpectrumTemplate(SpectrumType spType, ColorSpace space, RealType e0, RealType e1, RealType e2) {
            RealType x, y, brightness;
            switch (space) {
                case ColorSpace::sRGB_NonLinear: {
                    e0 = sRGB_degamma(e0);
                    e1 = sRGB_degamma(e1);
                    e2 = sRGB_degamma(e2);
                    // pass to the sRGB process.
                }
                case ColorSpace::sRGB: {
                    RealType RGB[3] = {e0, e1, e2};
                    RealType XYZ[3];
                    switch (spType) {
                        case SpectrumType::Reflectance:
                            sRGB_E_to_XYZ(RGB, XYZ);
                            break;
                        case SpectrumType::Illuminant:
                            sRGB_E_to_XYZ(RGB, XYZ);
                            break;
                        default:
                            SLRAssert(false, "Invalid Spectrum Type");
                            break;
                    }
                    e0 = XYZ[0];
                    e1 = XYZ[1];
                    e2 = XYZ[2];
                    // pass to the XYZ process.
                }
                case ColorSpace::XYZ: {
                    brightness = e0 + e1 + e2;
                    if (brightness == 0) {
                        u = 6;
                        v = 4;
                        scale = 0;
                        return;
                    }
                    x = e0 / brightness;
                    y = e1 / brightness;
                    break;
                }
                case ColorSpace::xyY: {
                    x = e0;
                    y = e1;
                    brightness = e2 / e1;
                    break;
                }
                default:
                    SLRAssert(false, "Invalid color space is specified.");
                    break;
            }
            // TODO: Contain a factor for solid of natural reflectance.
            //if (spType == SpectrumType::Reflectance)
            //    brightness = std::min(brightness, evaluateMaximumBrightness(x, y));
            scale = brightness / Upsampling::EqualEnergyReflectance;
            RealType xy[2] = {x, y};
            Upsampling::xy_to_uv(xy, &u);
            SLRAssert(!std::isinf(u) && !std::isnan(u) && !std::isinf(v) && !std::isnan(v) && !std::isinf(scale) && !std::isnan(scale), "Invalid value.");
        };
        
        SampledSpectrumTemplate<RealType, N> evaluate(const WavelengthSamplesTemplate<RealType, N> &wls) const override {
            using namespace Upsampling;
            if (u < 0.0f || u >= GridWidth || v < 0.0f || v >= GridHeight)
                return SampledSpectrumTemplate<RealType, N>::Zero;
            
            int32_t ui = (int32_t)u;
            int32_t vi = (int32_t)v;
            SLRAssert(ui < GridWidth && vi < GridHeight, "out of grid: %d, %d", ui, vi);
            
            const int32_t cellIdx = ui + GridWidth * vi;
            SLRAssert(cellIdx >= 0 && cellIdx < GridWidth * GridHeight, "cellIdx is out of grid: %d", cellIdx);
            
            const spectrum_grid_cell_t* cell = spectrum_grid + cellIdx;
            const uint8_t* indices = cell->idx;
            const uint8_t numPoints = cell->num_points;
            
            uint8_t usedIndices[4] = {UINT8_MAX, UINT8_MAX, UINT8_MAX, UINT8_MAX};
            float weights[4] = {0, 0, 0, 0};
            if (cell->inside) { // fast path for normal inner quads:
                // the layout of the vertices in the quad is:
                //  2  3
                //  0  1
                RealType s = u - ui;
                RealType t = v - vi;
                SLRAssert(s >= 0 && s <= 1 && t >= 0 && t <= 1, "invalid coordinate.");
                weights[0] = (1 - s) * (1 - t);
                weights[1] = s * (1 - t);
                weights[2] = (1 - s) * t;
                weights[3] = s * t;
                usedIndices[0] = indices[0];
                usedIndices[1] = indices[1];
                usedIndices[2] = indices[2];
                usedIndices[3] = indices[3];
            }
            else {
                // need to go through triangulation :(
                // we get the indices in such an order that they form a triangle fan around idx[0].
                // compute barycentric coordinates of our xy* point for all triangles in the fan:
                const float ex = u - spectrum_data_points[indices[0]].uv[0];
                const float ey = v - spectrum_data_points[indices[0]].uv[1];
                float e0x = spectrum_data_points[indices[1]].uv[0] - spectrum_data_points[indices[0]].uv[0];
                float e0y = spectrum_data_points[indices[1]].uv[1] - spectrum_data_points[indices[0]].uv[1];
                float uu = e0x * ey - ex * e0y;
                for (int i = 1; i < numPoints; ++i) {
                    uint32_t idx = indices[i % (numPoints - 1) + 1];
                    float e1x = spectrum_data_points[idx].uv[0] - spectrum_data_points[indices[0]].uv[0];
                    float e1y = spectrum_data_points[idx].uv[1] - spectrum_data_points[indices[0]].uv[1];
                    float vv = ex * e1y - e1x * ey;
                    
                    // TODO: with some sign magic, this division could be deferred to the last iteration!
                    const float area = e0x * e1y - e1x * e0y;
                    // normalise
                    const float u = uu / area;
                    const float v = vv / area;
                    float w = 1.0f - u - v;
                    // outside spectral locus (quantized version at least) or outside grid
                    if (u < 0.0 || v < 0.0 || w < 0.0) {
                        uu = -vv;
                        e0x = e1x;
                        e0y = e1y;
                        continue;
                    }
                    
                    weights[0] = u;
                    weights[1] = v;
                    weights[2] = w;
                    usedIndices[0] = idx;
                    usedIndices[1] = indices[i];
                    usedIndices[2] = indices[0];
                    break;
                }
            }
            
            SampledSpectrumTemplate<RealType, N> ret(0.0);
            for (int i = 0; i < WavelengthSamplesTemplate<RealType, N>::NumComponents; ++i) {
                RealType lambda = wls[i];
                RealType p = (lambda - MinWavelength) / (MaxWavelength - MinWavelength);
                SLRAssert(p >= 0 && p <= 1, "Wavelength is out of valid range.");
                RealType sBinF = p * (NumWavelengthSamples - 1);
                uint32_t sBin = (uint32_t)sBinF;
                uint32_t sBinNext = (sBin + 1 < NumWavelengthSamples) ? (sBin + 1) : (NumWavelengthSamples - 1);
                SLRAssert(sBin < NumWavelengthSamples && sBinNext < NumWavelengthSamples, "Spectrum bin index is out of range.");
                RealType t = sBinF - sBin;
                if (cell->inside) {
                    for (int j = 0; j < 4; ++j) {
                        const float* spectrum = spectrum_data_points[usedIndices[j]].spectrum;
                        ret[i] += weights[j] * (spectrum[sBin] * (1 - t) + spectrum[sBinNext] * t);
                    }
                }
                else {
                    for (int j = 0; j < 3; ++j) {
                        const float* spectrum = spectrum_data_points[usedIndices[j]].spectrum;
                        ret[i] += weights[j] * (spectrum[sBin] * (1 - t) + spectrum[sBinNext] * t);
                    }
                }
            }
            
            return ret * scale;
        };
    };


    template <typename RealType, uint32_t N>
    struct SLR_API SampledSpectrumTemplate {
        RealType values[N];
        static const uint32_t NumComponents;
        
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
            RealType r = 1 / s;
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
            SLRAssert(index < N, "\"index\" is out of range [0, %u].", N - 1);
            return values[index];
        };
        RealType operator[](unsigned int index) const {
            SLRAssert(index < N, "\"index\" is out of range [0, %u].", N - 1);
            return values[index];
        };
        
        RealType avgValue() const {
            RealType sumVal = values[0];
            for (int i = 1; i < N; ++i)
                sumVal += values[i];
            return sumVal / N;
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
        bool hasMinus() const {
            for (int i = 0; i < N; ++i)
                if (values[i] < 0)
                    return true;
            return false;
        };
        
        RealType luminance(RGBColorSpace space = RGBColorSpace::sRGB) const {
            RealType sum = 0;
            for (int i = 0; i < N; ++i)
                sum += values[i];
            return sum;
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
            return ret;
        };
        
        std::string toString(const WavelengthSamplesTemplate<RealType, N> &wls) const {
            std::string ret = "";
            char str[256];
            for (int i = 0; i < N; ++i) {
                sprintf(str, "%g, %g\n", wls[i], values[i]);
                ret += str;
            }
            return ret;
        };
        
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


    template <typename RealType, uint32_t numStrata>
    struct SLR_API DiscretizedSpectrumTemplate {
        RealType values[numStrata];
        
        DiscretizedSpectrumTemplate(RealType v = 0.0f) { for (int i = 0; i < numStrata; ++i) values[i] = v; };
        DiscretizedSpectrumTemplate(const RealType* vals) { for (int i = 0; i < numStrata; ++i) values[i] = vals[i]; };
        
        DiscretizedSpectrumTemplate operator+() const { return *this; };
        DiscretizedSpectrumTemplate operator-() const {
            RealType vals[numStrata];
            for (int i = 0; i < numStrata; ++i)
                vals[i] = -values[i];
            return DiscretizedSpectrumTemplate(vals);
        };
        
        DiscretizedSpectrumTemplate operator+(const DiscretizedSpectrumTemplate &c) const {
            RealType vals[numStrata];
            for (int i = 0; i < numStrata; ++i)
                vals[i] = values[i] + c.values[i];
            return DiscretizedSpectrumTemplate(vals);
        };
        DiscretizedSpectrumTemplate operator-(const DiscretizedSpectrumTemplate &c) const {
            RealType vals[numStrata];
            for (int i = 0; i < numStrata; ++i)
                vals[i] = values[i] - c.values[i];
            return DiscretizedSpectrumTemplate(vals);
        };
        DiscretizedSpectrumTemplate operator*(const DiscretizedSpectrumTemplate &c) const {
            RealType vals[numStrata];
            for (int i = 0; i < numStrata; ++i)
                vals[i] = values[i] * c.values[i];
            return DiscretizedSpectrumTemplate(vals);
        };
        DiscretizedSpectrumTemplate operator*(RealType s) const {
            RealType vals[numStrata];
            for (int i = 0; i < numStrata; ++i)
                vals[i] = values[i] * s;
            return DiscretizedSpectrumTemplate(vals);
        };
        friend inline DiscretizedSpectrumTemplate operator*(RealType s, const DiscretizedSpectrumTemplate &c) {
            RealType vals[numStrata];
            for (int i = 0; i < numStrata; ++i)
                vals[i] = c.values[i] * s;
            return DiscretizedSpectrumTemplate(vals);
        };
        
        DiscretizedSpectrumTemplate &operator+=(const DiscretizedSpectrumTemplate &c) {
            for (int i = 0; i < numStrata; ++i)
                values[i] += c.values[i];
            return *this;
        };
        DiscretizedSpectrumTemplate &operator*=(const DiscretizedSpectrumTemplate &c) {
            for (int i = 0; i < numStrata; ++i)
                values[i] *= c.values[i];
            return *this;
        };
        DiscretizedSpectrumTemplate &operator*=(RealType s) {
            for (int i = 0; i < numStrata; ++i)
                values[i] *= s;
            return *this;
        };
        
        bool operator==(const DiscretizedSpectrumTemplate &c) const {
            for (int i = 0; i < numStrata; ++i)
                if (values[i] != c.values[i])
                    return false;
            return true;
        };
        bool operator!=(const DiscretizedSpectrumTemplate &c) const {
            for (int i = 0; i < numStrata; ++i)
                if (values[i] != c.values[i])
                    return true;
            return false;
        };
        
        RealType &operator[](unsigned int index) {
            SLRAssert(index < numStrata, "\"index\" is out of range [0, %u].", numStrata - 1);
            return values[index];
        };
        RealType operator[](unsigned int index) const {
            SLRAssert(index < numStrata, "\"index\" is out of range [0, %u].", numStrata - 1);
            return values[index];
        };
        
        RealType maxValue() const {
            RealType maxVal = values[0];
            for (int i = 1; i < numStrata; ++i)
                maxVal = std::fmax(values[i], maxVal);
            return maxVal;
        };
        RealType minValue() const {
            RealType minVal = values[0];
            for (int i = 1; i < numStrata; ++i)
                minVal = std::fmin(values[i], minVal);
            return minVal;
        };
        bool hasNonZero() const {
            for (int i = 0; i < numStrata; ++i)
                if (values[i] != 0)
                    return true;
            return false;
        };
        bool hasNaN() const {
            for (int i = 0; i < numStrata; ++i)
                if (std::isnan(values[i]))
                    return true;
            return false;
        };
        bool hasInf() const {
            for (int i = 0; i < numStrata; ++i)
                if (std::isinf(values[i]))
                    return true;
            return false;
        };
        bool hasMinus() const {
            for (int i = 0; i < numStrata; ++i)
                if (values[i] < 0)
                    return true;
            return false;
        };
        
        RealType luminance(RGBColorSpace space = RGBColorSpace::sRGB) const {
            RealType sum = 0;
            for (int i = 0; i < numStrata; ++i)
                sum += ybar[i] * values[i];
            return sum / integralCMF;
        };
        
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
        };
        
        std::string toString() const {
            std::string ret = "(";
            char str[256];
            for (int i = 0; i < numStrata - 1; ++i) {
                sprintf(str, "%g, ", values[i]);
                ret += str;
            }
            sprintf(str, "%g)", values[numStrata - 1]);
            ret += str;
            return str;
        };
        
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
        };
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
    struct SLR_API SpectrumStorageTemplate {
        typedef DiscretizedSpectrumTemplate<RealType, numStrata> ValueType;
        CompensatedSum<ValueType> value;
        
        SpectrumStorageTemplate(const ValueType &v = ValueType::Zero) :
        value(v) {};
        
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
        };
    };
}

#endif
