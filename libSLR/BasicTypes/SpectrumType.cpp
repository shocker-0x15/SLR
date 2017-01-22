//
//  SpectrumType.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "SpectrumTypes.h"

namespace SLR {
    template struct SLR_API WavelengthSamplesTemplate<float, NumSpectralSamples>;
//    template SLR_API const uint32_t WavelengthSamplesTemplate<float, NumSpectralSamples>::NumComponents;

    template struct SLR_API WavelengthSamplesTemplate<double, NumSpectralSamples>;
//    template SLR_API const uint32_t WavelengthSamplesTemplate<double, NumSpectralSamples>::NumComponents;

    
    template struct SLR_API ContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct SLR_API ContinuousSpectrumTemplate<double, NumSpectralSamples>;


    
    template <typename RealType, uint32_t N>
    RealType RegularContinuousSpectrumTemplate<RealType, N>::calcBounds() const {
        RealType maxValue = -INFINITY;
        for (int i = 0; i < numSamples; ++i)
            maxValue = std::max(maxValue, values[i]);
        return maxValue;
    }
    
    template <typename RealType, uint32_t N>
    SampledSpectrumTemplate<RealType, N> RegularContinuousSpectrumTemplate<RealType, N>::evaluate(const WavelengthSamplesTemplate<RealType, N> &wls) const {
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
    }
    
    template <typename RealType, uint32_t N>
    ContinuousSpectrumTemplate<RealType, N>* RegularContinuousSpectrumTemplate<RealType, N>::createScaled(RealType scale) const {
        RealType* sValues = new RealType[numSamples];
        for (int i = 0; i < numSamples; ++i)
            sValues[i] = scale * values[i];
        ContinuousSpectrumTemplate<RealType, N>* ret = new RegularContinuousSpectrumTemplate(minLambda, maxLambda, sValues, numSamples);
        delete[] sValues;
        return ret;
    }
    
    template struct SLR_API RegularContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct SLR_API RegularContinuousSpectrumTemplate<double, NumSpectralSamples>;

    
    
    template <typename RealType, uint32_t N>
    RealType IrregularContinuousSpectrumTemplate<RealType, N>::calcBounds() const {
        RealType maxValue = -INFINITY;
        for (int i = 0; i < numSamples; ++i)
            maxValue = std::max(maxValue, values[i]);
        return maxValue;
    }
    
    template <typename RealType, uint32_t N>
    SampledSpectrumTemplate<RealType, N> IrregularContinuousSpectrumTemplate<RealType, N>::evaluate(const WavelengthSamplesTemplate<RealType, N> &wls) const {
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
    }
    
    template <typename RealType, uint32_t N>
    ContinuousSpectrumTemplate<RealType, N>* IrregularContinuousSpectrumTemplate<RealType, N>::createScaled(RealType scale) const {
        RealType* sValues = new RealType[numSamples];
        for (int i = 0; i < numSamples; ++i)
            sValues[i] = scale * values[i];
        ContinuousSpectrumTemplate<RealType, N>* ret = new IrregularContinuousSpectrumTemplate(lambdas, sValues, numSamples);
        delete[] sValues;
        return ret;
    }
    
    template struct SLR_API IrregularContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct SLR_API IrregularContinuousSpectrumTemplate<double, NumSpectralSamples>;

    
    
    template <typename RealType, uint32_t N>
    void UpsampledContinuousSpectrumTemplate<RealType, N>::computeAdjacents(RealType u, RealType v) {
        using namespace Upsampling;
        u = std::clamp<RealType>(u, 0.0f, GridWidth);
        v = std::clamp<RealType>(v, 0.0f, GridHeight);
        
        int32_t ui = (int32_t)u;
        int32_t vi = (int32_t)v;
        SLRAssert(ui < GridWidth && vi < GridHeight, "out of grid: %d, %d", ui, vi);
        
        const int32_t cellIdx = ui + GridWidth * vi;
        SLRAssert(cellIdx >= 0 && cellIdx < GridWidth * GridHeight, "cellIdx is out of grid: %d", cellIdx);
        
        const spectrum_grid_cell_t* cell = spectrum_grid + cellIdx;
        const uint8_t* indices = cell->idx;
        const uint8_t numPoints = cell->num_points;
        
        if (cell->inside) { // fast path for normal inner quads:
            // the layout of the vertices in the quad is:
            //  2  3
            //  0  1
            RealType s = u - ui;
            RealType t = v - vi;
            SLRAssert(s >= 0 && s <= 1 && t >= 0 && t <= 1, "invalid coordinate.");
            m_adjIndices = (((uint32_t)indices[3] << 24) |
                            ((uint32_t)indices[2] << 16) |
                            ((uint32_t)indices[1] << 8) |
                            ((uint32_t)indices[0] << 0));
            m_s = s;
            m_t = t;
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
                const float b0 = uu / area;
                const float b1 = vv / area;
                float b2 = 1.0f - b0 - b1;
                // outside spectral locus (quantized version at least) or outside grid
                //if (b0 < 0.0 || b1 < 0.0 || b2 < 0.0) {
                if (b0 < -1e-6 || b1 < -1e-6 || b2 < -1e-6) {
                    uu = -vv;
                    e0x = e1x;
                    e0y = e1y;
                    continue;
                }
                
                m_adjIndices = (((uint32_t)UINT8_MAX << 24) |
                                ((uint32_t)indices[0] << 16) |
                                ((uint32_t)indices[i] << 8) |
                                ((uint32_t)idx << 0));
                m_s = b0;
                m_t = b1;
                break;
            }
        }
        SLRAssert((m_adjIndices && 0xFF) != UINT8_MAX, "Adjacent points must be selected at this point.");
    }
    
    template <typename RealType, uint32_t N>
    UpsampledContinuousSpectrumTemplate<RealType, N>::UpsampledContinuousSpectrumTemplate(SpectrumType spType, ColorSpace space, RealType e0, RealType e1, RealType e2) {
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
                        sRGB_to_XYZ(RGB, XYZ);
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
                    m_scale = 0;
                    computeAdjacents(6, 4);
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
        m_scale = brightness / Upsampling::EqualEnergyReflectance;
        RealType xy[2] = {x, y};
        RealType uv[2];
        Upsampling::xy_to_uv(xy, uv);
        SLRAssert(!std::isinf(uv[0]) && !std::isnan(uv[0]) &&
                  !std::isinf(uv[1]) && !std::isnan(uv[1]) &&
                  !std::isinf(m_scale) && !std::isnan(m_scale), "Invalid value.");
        
        computeAdjacents(uv[0], uv[1]);
    }
    
    template <typename RealType, uint32_t N>
    RealType UpsampledContinuousSpectrumTemplate<RealType, N>::calcBounds() const {
        using namespace Upsampling;
        
        uint8_t adjIndices[4];
        adjIndices[0] = (m_adjIndices >> 0) & 0xFF;
        adjIndices[1] = (m_adjIndices >> 8) & 0xFF;
        adjIndices[2] = (m_adjIndices >> 16) & 0xFF;
        adjIndices[3] = (m_adjIndices >> 24) & 0xFF;
        
        int numAdjacents = (adjIndices[3] != UINT8_MAX) ? 4 : 3;
        
        RealType maxValue = -INFINITY;
        for (int i = 0; i < numAdjacents; ++i) {
            const float* spectrum = spectrum_data_points[adjIndices[i]].spectrum;
            for (int j = 0; j < NumWavelengthSamples; ++j)
                maxValue = std::max<RealType>(maxValue, spectrum[j]);
        }
        
        return maxValue;
    }
    
    template <typename RealType, uint32_t N>
    SampledSpectrumTemplate<RealType, N> UpsampledContinuousSpectrumTemplate<RealType, N>::evaluate(const WavelengthSamplesTemplate<RealType, N> &wls) const {
        using namespace Upsampling;
        
        uint8_t adjIndices[4];
        adjIndices[0] = (m_adjIndices >> 0) & 0xFF;
        adjIndices[1] = (m_adjIndices >> 8) & 0xFF;
        adjIndices[2] = (m_adjIndices >> 16) & 0xFF;
        adjIndices[3] = (m_adjIndices >> 24) & 0xFF;
        
        int numAdjacents;
        float weights[4];
        if (adjIndices[3] != UINT8_MAX) {
            weights[0] = (1 - m_s) * (1 - m_t);
            weights[1] = m_s * (1 - m_t);
            weights[2] = (1 - m_s) * m_t;
            weights[3] = m_s * m_t;
            numAdjacents = 4;
        }
        else {
            weights[0] = m_s;
            weights[1] = m_t;
            weights[2] = 1.0f - m_s - m_t;
            numAdjacents = 3;
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
            for (int j = 0; j < numAdjacents; ++j) {
                const float* spectrum = spectrum_data_points[adjIndices[j]].spectrum;
                ret[i] += weights[j] * (spectrum[sBin] * (1 - t) + spectrum[sBinNext] * t);
            }
        }
        
        return ret * m_scale;
    }
    
    template struct SLR_API UpsampledContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct SLR_API UpsampledContinuousSpectrumTemplate<double, NumSpectralSamples>;
    

    template struct SLR_API SampledSpectrumTemplate<float, NumSpectralSamples>;
//    template SLR_API const uint32_t SampledSpectrumTemplate<float, NumSpectralSamples>::NumComponents;
//    template SLR_API const SampledSpectrumTemplate<float, NumSpectralSamples> SampledSpectrumTemplate<float, NumSpectralSamples>::Zero;
//    template SLR_API const SampledSpectrumTemplate<float, NumSpectralSamples> SampledSpectrumTemplate<float, NumSpectralSamples>::One;
//    template SLR_API const SampledSpectrumTemplate<float, NumSpectralSamples> SampledSpectrumTemplate<float, NumSpectralSamples>::Inf;
//    template SLR_API const SampledSpectrumTemplate<float, NumSpectralSamples> SampledSpectrumTemplate<float, NumSpectralSamples>::NaN;

    template struct SLR_API SampledSpectrumTemplate<double, NumSpectralSamples>;
//    template SLR_API const uint32_t SampledSpectrumTemplate<double, NumSpectralSamples>::NumComponents;
//    template SLR_API const SampledSpectrumTemplate<double, NumSpectralSamples> SampledSpectrumTemplate<double, NumSpectralSamples>::Zero;
//    template SLR_API const SampledSpectrumTemplate<double, NumSpectralSamples> SampledSpectrumTemplate<double, NumSpectralSamples>::One;
//    template SLR_API const SampledSpectrumTemplate<double, NumSpectralSamples> SampledSpectrumTemplate<double, NumSpectralSamples>::Inf;
//    template SLR_API const SampledSpectrumTemplate<double, NumSpectralSamples> SampledSpectrumTemplate<double, NumSpectralSamples>::NaN;
    
    template <typename RealType, uint32_t N>
    SampledSpectrumTemplate<RealType, N> sqrt(const SampledSpectrumTemplate<RealType, N> &value) {
        SampledSpectrumTemplate<RealType, N> ret;
        for (int i = 0; i < N; ++i)
            ret[i] = std::sqrt(value[i]);
        return ret;
    }
    template SLR_API SampledSpectrumTemplate<float, NumSpectralSamples> sqrt(const SampledSpectrumTemplate<float, NumSpectralSamples> &value);
    template SLR_API SampledSpectrumTemplate<double, NumSpectralSamples> sqrt(const SampledSpectrumTemplate<double, NumSpectralSamples> &value);
    
    template <typename RealType, uint32_t N>
    SampledSpectrumTemplate<RealType, N> exp(const SampledSpectrumTemplate<RealType, N> &value) {
        SampledSpectrumTemplate<RealType, N> ret;
        for (int i = 0; i < N; ++i)
            ret[i] = std::exp(value[i]);
        return ret;
    }
    template SLR_API SampledSpectrumTemplate<float, NumSpectralSamples> exp(const SampledSpectrumTemplate<float, NumSpectralSamples> &value);
    template SLR_API SampledSpectrumTemplate<double, NumSpectralSamples> exp(const SampledSpectrumTemplate<double, NumSpectralSamples> &value);

    
    template struct DiscretizedSpectrumTemplate<float, NumStrataForStorage>;
//    template SLR_API const uint32_t DiscretizedSpectrumTemplate<float, NumStrataForStorage>::NumStrata;
//    template SLR_API std::unique_ptr<float[]> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::xbar;
//    template SLR_API std::unique_ptr<float[]> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::ybar;
//    template SLR_API std::unique_ptr<float[]> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::zbar;
//    template SLR_API float DiscretizedSpectrumTemplate<float, NumStrataForStorage>::integralCMF;
//    template SLR_API const DiscretizedSpectrumTemplate<float, NumStrataForStorage> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::Zero;
//    template SLR_API const DiscretizedSpectrumTemplate<float, NumStrataForStorage> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::One;
//    template SLR_API const DiscretizedSpectrumTemplate<float, NumStrataForStorage> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::Inf;
//    template SLR_API const DiscretizedSpectrumTemplate<float, NumStrataForStorage> DiscretizedSpectrumTemplate<float, NumStrataForStorage>::NaN;

    template struct DiscretizedSpectrumTemplate<double, NumStrataForStorage>;
//    template SLR_API const uint32_t DiscretizedSpectrumTemplate<double, NumStrataForStorage>::NumStrata;
//    template SLR_API std::unique_ptr<double[]> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::xbar;
//    template SLR_API std::unique_ptr<double[]> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::ybar;
//    template SLR_API std::unique_ptr<double[]> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::zbar;
//    template SLR_API double DiscretizedSpectrumTemplate<double, NumStrataForStorage>::integralCMF;
//    template SLR_API const DiscretizedSpectrumTemplate<double, NumStrataForStorage> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::Zero;
//    template SLR_API const DiscretizedSpectrumTemplate<double, NumStrataForStorage> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::One;
//    template SLR_API const DiscretizedSpectrumTemplate<double, NumStrataForStorage> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::Inf;
//    template SLR_API const DiscretizedSpectrumTemplate<double, NumStrataForStorage> DiscretizedSpectrumTemplate<double, NumStrataForStorage>::NaN;
    
    template struct SLR_API SpectrumStorageTemplate<float, NumStrataForStorage>;
    template struct SLR_API SpectrumStorageTemplate<double, NumStrataForStorage>;
}
