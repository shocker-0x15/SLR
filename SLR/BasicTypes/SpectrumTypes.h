//
//  SpectrumTypes.h
//
//  Created by 渡部 心 on 2015/09/13.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef SLR_SpectrumTypes_h
#define SLR_SpectrumTypes_h

template <typename RealType, uint32_t N>
struct WavelengthSamplesTemplate {
    RealType lambdas[N];
    uint32_t flags;
    static const uint32_t NumComponents;

    RealType &operator[](uint32_t index) {
        SLRAssert(index < N, "\"index\" is out of range [0, %u].", N - 1);
        return lambdas[index];
    };
    RealType operator[](uint32_t index) const {
        SLRAssert(index < N, "\"index\" is out of range [0, %u].", N - 1);
        return lambdas[index];
    };
};
template <typename RealType, uint32_t N> const uint32_t WavelengthSamplesTemplate<RealType, N>::NumComponents = N;

template <typename RealType, uint32_t N>
struct ContinuousSpectrumTemplate {
    virtual SampledSpectrumTemplate<RealType, N> evaluate(const WavelengthSamples &wls) const = 0;
};

template <typename RealType, uint32_t N>
struct RegularContinuousSpectrumTemplate : public ContinuousSpectrumTemplate<RealType, N> {
    SampledSpectrumTemplate<RealType, N> evaluate(const WavelengthSamples &wls) const override {
    
    };
};

template <typename RealType, uint32_t N>
struct UpsampledContinuousSpectrumTemplate : public ContinuousSpectrumTemplate<RealType, N> {
    RealType u, v, norm;
    
    std::function<RealType(RealType)> inverseGammaCorrection = [](RealType value) {
        if (value <= 0.04045)
            return value  / 12.92;
        return std::pow((value + 0.055) / 1.055, 2.4);
    };
    
    UpsampledContinuousSpectrumTemplate(ColorSpace space, RealType e0, RealType e1, RealType e2) {
        RealType x, y;
        switch (space) {
            case ColorSpace::sRGB_NonLinear: {
                e0 = inverseGammaCorrection(e0);
                e1 = inverseGammaCorrection(e1);
                e2 = inverseGammaCorrection(e2);
                // pass to the sRGB process.
            }
            case ColorSpace::sRGB: {
                RealType RGB[3] = {e0, e1, e2};
                RealType XYZ[3];
                sRGB_to_XYZ(RGB, XYZ);
                e0 = XYZ[0];
                e1 = XYZ[1];
                e2 = XYZ[2];
                // pass to the XYZ process.
            }
            case ColorSpace::XYZ: {
                norm = e0 + e1 + e2;
                x = e0 / norm;
                y = e1 / norm;
                break;
            }
            case ColorSpace::xyY: {
                x = e0;
                y = e1;
                norm = e2 / e1;
                break;
            }
            default:
                SLRAssert(false, "Invalid color space is specified.");
                break;
        }
        RealType xy[2] = {x, y};
        xy_to_uv(xy, &u);
    };
    
    SampledSpectrumTemplate<RealType, N> evaluate(const WavelengthSamples &wls) const override {
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
        for (int i = 0; i < WavelengthSamples::NumComponents; ++i) {
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

        return ret * norm / EqualEnergyReflectance;
    };
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
        SLRAssert(index < N, "\"index\" is out of range [0, %u].", N - 1);
        return values[index];
    };
    RealType operator[](unsigned int index) const {
        SLRAssert(index < N, "\"index\" is out of range [0, %u].", N - 1);
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
    
    RealType luminance(RGBColorSpace space = RGBColorSpace::sRGB) const {
        RealType sum = 0;
        for (int i = 0; i < N; ++i)
            sum += i;
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
template <typename RealType, uint32_t N> const uint32_t SampledSpectrumTemplate<RealType, N>::numComponents = N;

template <typename RealType, uint32_t N>
const SampledSpectrumTemplate<RealType, N> SampledSpectrumTemplate<RealType, N>::Zero = SampledSpectrumTemplate<RealType, N>(0.0);
template <typename RealType, uint32_t N>
const SampledSpectrumTemplate<RealType, N> SampledSpectrumTemplate<RealType, N>::One = SampledSpectrumTemplate<RealType, N>(1.0);
template <typename RealType, uint32_t N>
const SampledSpectrumTemplate<RealType, N> SampledSpectrumTemplate<RealType, N>::Inf = SampledSpectrumTemplate<RealType, N>(std::numeric_limits<RealType>::infinity());
template <typename RealType, uint32_t N>
const SampledSpectrumTemplate<RealType, N> SampledSpectrumTemplate<RealType, N>::NaN = SampledSpectrumTemplate<RealType, N>(std::numeric_limits<RealType>::quiet_NaN());


template <typename RealType, uint32_t numStrata>
struct DiscretizedSpectrumTemplate {
    RealType values[numStrata];
    static const uint32_t NumStrata;
    
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
    
    float luminance(RGBColorSpace space = RGBColorSpace::sRGB) const {
        return 0.0f;
    };
    
    void getRGB(float RGB[3], RGBColorSpace space = RGBColorSpace::sRGB) const {
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
};
template <typename RealType, uint32_t numStrata> const uint32_t DiscretizedSpectrumTemplate<RealType, numStrata>::NumStrata = numStrata;

template <typename RealType, uint32_t numStrata>
const DiscretizedSpectrumTemplate<RealType, numStrata> DiscretizedSpectrumTemplate<RealType, numStrata>::Zero = DiscretizedSpectrumTemplate<RealType, numStrata>(0.0);
template <typename RealType, uint32_t numStrata>
const DiscretizedSpectrumTemplate<RealType, numStrata> DiscretizedSpectrumTemplate<RealType, numStrata>::One = DiscretizedSpectrumTemplate<RealType, numStrata>(1.0);
template <typename RealType, uint32_t numStrata>
const DiscretizedSpectrumTemplate<RealType, numStrata> DiscretizedSpectrumTemplate<RealType, numStrata>::Inf = DiscretizedSpectrumTemplate<RealType, numStrata>(std::numeric_limits<RealType>::infinity());
template <typename RealType, uint32_t numStrata>
const DiscretizedSpectrumTemplate<RealType, numStrata> DiscretizedSpectrumTemplate<RealType, numStrata>::NaN = DiscretizedSpectrumTemplate<RealType, numStrata>(std::numeric_limits<RealType>::quiet_NaN());


template <typename RealType, uint32_t numStrata>
struct SpectrumStorageTemplate {
    CompensatedSum<DiscretizedSpectrumTemplate<RealType, numStrata>> value;
    
    SpectrumStorageTemplate(const DiscretizedSpectrumTemplate<RealType, numStrata> &v = DiscretizedSpectrumTemplate<RealType, numStrata>::Zero) :
    value(v) {};
    
    template <uint32_t N>
    SpectrumStorageTemplate &add(const WavelengthSamplesTemplate<RealType, N> &wls, const SampledSpectrumTemplate<RealType, N> &value) {
        return *this;
    };
};

#endif
