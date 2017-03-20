//
//  spectrum_types.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "spectrum_types.h"

namespace SLR {
    template struct SLR_API WavelengthSamplesTemplate<float, NumSpectralSamples>;
//    template SLR_API const uint32_t WavelengthSamplesTemplate<float, NumSpectralSamples>::NumComponents;

    template struct SLR_API WavelengthSamplesTemplate<double, NumSpectralSamples>;
//    template SLR_API const uint32_t WavelengthSamplesTemplate<double, NumSpectralSamples>::NumComponents;

    
    template class SLR_API ContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template class SLR_API ContinuousSpectrumTemplate<double, NumSpectralSamples>;


    
    template <typename RealType, uint32_t NumSpectralSamples>
    void RegularContinuousSpectrumTemplate<RealType, NumSpectralSamples>::calcBounds(uint32_t numBins, RealType* bounds) const {
        const RealType BinWidth = (WavelengthHighBound - WavelengthLowBound) / numBins;
        const RealType SampleInterval = (m_maxLambda - m_minLambda) / (m_numSamples - 1);
        
        for (int binIdx = 0; binIdx < numBins; ++binIdx) {
            RealType wlLowInBin = WavelengthLowBound + BinWidth * binIdx;
            RealType wlHighInBin = WavelengthLowBound + BinWidth * (binIdx + 1);
            
            int lowestWlIdxInBin = std::max(0, (int)((wlLowInBin - m_minLambda) / SampleInterval) + 1);
            int highestWlIdxInBin = std::min((int)(m_numSamples - 1), (int)((wlHighInBin - m_minLambda) / SampleInterval));
            int highestWlIdxOutBin = std::max(0, lowestWlIdxInBin - 1);
            int lowestWlIdxOutBin = std::min((int)(m_numSamples - 1), highestWlIdxInBin + 1);
            RealType lerpParamAtLowBoundary = (wlLowInBin - (m_minLambda + highestWlIdxOutBin * SampleInterval)) / SampleInterval;
            RealType lerpParamAtHighBoundary = (wlHighInBin - (m_minLambda + highestWlIdxInBin * SampleInterval)) / SampleInterval;
            
            RealType maxValue = -INFINITY;
            RealType lerpedValueAtLow = m_values[highestWlIdxOutBin] * (1 - lerpParamAtLowBoundary) + m_values[lowestWlIdxInBin] * lerpParamAtLowBoundary;
            maxValue = std::max(maxValue, lerpedValueAtLow);
            for (int wlIdx = lowestWlIdxInBin; wlIdx <= highestWlIdxInBin; ++wlIdx)
                maxValue = std::max(maxValue, m_values[wlIdx]);
            RealType lerpedValueAtHigh = m_values[highestWlIdxInBin] * (1 - lerpParamAtHighBoundary) + m_values[lowestWlIdxOutBin] * lerpParamAtHighBoundary;
            maxValue = std::max(maxValue, lerpedValueAtHigh);
            
            bounds[binIdx] = maxValue;
        }
    }
    
    template <typename RealType, uint32_t NumSpectralSamples>
    SampledSpectrumTemplate<RealType, NumSpectralSamples> 
    RegularContinuousSpectrumTemplate<RealType, NumSpectralSamples>::evaluate(const WavelengthSamplesTemplate<RealType, NumSpectralSamples> &wls) const {
        SampledSpectrumTemplate<RealType, NumSpectralSamples> ret(0.0f);
        for (int i = 0; i < WavelengthSamplesTemplate<RealType, NumSpectralSamples>::NumComponents; ++i) {
            RealType binF = (wls[i] - m_minLambda) / (m_maxLambda - m_minLambda) * (m_numSamples - 1);
            if (binF <= 0.0f) {
                ret[i] = m_values[0];
                continue;
            }
            else if (binF >= m_numSamples - 1) {
                ret[i] = m_values[m_numSamples - 1];
                continue;
            }
            int32_t bin = int32_t(binF);
            SLRAssert(bin >= 0 && bin < m_numSamples - 1, "invalid bin index.");
            RealType t = binF - bin;
            ret[i] = (1 - t) * m_values[bin] + t * m_values[bin + 1];
        }
        return ret;
    }
    
    template <typename RealType, uint32_t NumSpectralSamples>
    ContinuousSpectrumTemplate<RealType, NumSpectralSamples>* 
    RegularContinuousSpectrumTemplate<RealType, NumSpectralSamples>::createScaledAndOffset(RealType scale, RealType offset) const {
        RealType* sValues = new RealType[m_numSamples];
        for (int i = 0; i < m_numSamples; ++i)
            sValues[i] = scale * m_values[i] + offset;
        ContinuousSpectrumTemplate<RealType, NumSpectralSamples>* ret = new RegularContinuousSpectrumTemplate(m_minLambda, m_maxLambda, sValues, m_numSamples);
        delete[] sValues;
        return ret;
    }
    
    template class SLR_API RegularContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template class SLR_API RegularContinuousSpectrumTemplate<double, NumSpectralSamples>;

    
    
    template <typename RealType, uint32_t NumSpectralSamples>
    void IrregularContinuousSpectrumTemplate<RealType, NumSpectralSamples>::calcBounds(uint32_t numBins, RealType* bounds) const {
        const RealType BinWidth = (WavelengthHighBound - WavelengthLowBound) / numBins;
        
        int lastWlIdx = 0;
        for (int binIdx = 0; binIdx < numBins; ++binIdx) {
            RealType wlLowInBin = WavelengthLowBound + BinWidth * binIdx;
            RealType wlHighInBin = WavelengthLowBound + BinWidth * (binIdx + 1);
            
            int highestWlIdxOutBin = 0;
            for (int wlIdx = lastWlIdx; wlIdx < m_numSamples - 1; ++wlIdx) {
                if (m_lambdas[wlIdx + 1] >= wlLowInBin) {
                    highestWlIdxOutBin = wlIdx;
                    lastWlIdx = wlIdx;
                    break;
                }
            }
            RealType highestWlOutBin = m_lambdas[highestWlIdxOutBin];
            int lowestWlIdxInBin = std::min((int)m_numSamples - 1, highestWlIdxOutBin + 1);
            RealType lowestWlInBin = m_lambdas[lowestWlIdxInBin];
            
            int highestWlIdxInBin = 0;
            for (int wlIdx = lastWlIdx; wlIdx < m_numSamples - 1; ++wlIdx) {
                if (m_lambdas[wlIdx + 1] >= wlHighInBin) {
                    highestWlIdxInBin = wlIdx;
                    lastWlIdx = wlIdx;
                    break;
                }
            }
            RealType highestWlInBin = m_lambdas[highestWlIdxInBin];
            int lowestWlIdxOutBin = std::min((int)(m_numSamples - 1), highestWlIdxInBin + 1);
            RealType lowestWlOutBin = m_lambdas[lowestWlIdxOutBin]; 
            
            RealType lerpParamAtLowBoundary = (wlLowInBin - highestWlOutBin) / (lowestWlInBin - highestWlOutBin);
            RealType lerpParamAtHighBoundary = (wlHighInBin - highestWlInBin) / (lowestWlOutBin - highestWlInBin);
            if (lowestWlIdxInBin == highestWlIdxOutBin)
                lerpParamAtLowBoundary = 0.0f;
            if (highestWlIdxInBin == lowestWlIdxOutBin)
                lerpParamAtHighBoundary = 0.0f;
            
            RealType maxValue = -INFINITY;
            RealType lerpedValueAtLow = m_values[highestWlIdxOutBin] * (1 - lerpParamAtLowBoundary) + m_values[lowestWlIdxInBin] * lerpParamAtLowBoundary;
            maxValue = std::max(maxValue, lerpedValueAtLow);
            for (int wlIdx = lowestWlIdxInBin; wlIdx <= highestWlIdxInBin; ++wlIdx)
                maxValue = std::max(maxValue, m_values[wlIdx]);
            RealType lerpedValueAtHigh = m_values[highestWlIdxInBin] * (1 - lerpParamAtHighBoundary) + m_values[lowestWlIdxOutBin] * lerpParamAtHighBoundary;
            maxValue = std::max(maxValue, lerpedValueAtHigh);
            
            bounds[binIdx] = maxValue;
        }
    }
    
    template <typename RealType, uint32_t NumSpectralSamples>
    SampledSpectrumTemplate<RealType, NumSpectralSamples> 
    IrregularContinuousSpectrumTemplate<RealType, NumSpectralSamples>::evaluate(const WavelengthSamplesTemplate<RealType, NumSpectralSamples> &wls) const {
        SampledSpectrumTemplate<RealType, NumSpectralSamples> ret(0.0f);
        uint32_t searchBase = 0;
        for (int i = 0; i < WavelengthSamplesTemplate<RealType, NumSpectralSamples>::NumComponents; ++i) {
            int32_t lowIdx = std::max((int32_t)std::distance(m_lambdas, std::lower_bound(m_lambdas + searchBase, m_lambdas + m_numSamples, wls[i])) - 1, 0);
            searchBase = lowIdx;
            if (lowIdx >= m_numSamples - 1) {
                ret[i] = m_values[m_numSamples - 1];
                continue;
            }
            RealType t = (wls[i] - m_lambdas[lowIdx]) / (m_lambdas[lowIdx + 1] - m_lambdas[lowIdx]);
            if (t <= 0.0f) {
                ret[i] = m_values[0];
                continue;
            }
            SLRAssert(t >= 0 && t <= 1, "invalid interpolation coefficient.");
            ret[i] = (1 - t) * m_values[lowIdx] + t * m_values[lowIdx + 1];
        }
        return ret;
    }
    
    template <typename RealType, uint32_t NumSpectralSamples>
    ContinuousSpectrumTemplate<RealType, NumSpectralSamples>* 
    IrregularContinuousSpectrumTemplate<RealType, NumSpectralSamples>::createScaledAndOffset(RealType scale, RealType offset) const {
        RealType* sValues = new RealType[m_numSamples];
        for (int i = 0; i < m_numSamples; ++i)
            sValues[i] = scale * m_values[i] + offset;
        ContinuousSpectrumTemplate<RealType, NumSpectralSamples>* ret = new IrregularContinuousSpectrumTemplate(m_lambdas, sValues, m_numSamples);
        delete[] sValues;
        return ret;
    }
    
    template class SLR_API IrregularContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template class SLR_API IrregularContinuousSpectrumTemplate<double, NumSpectralSamples>;

    
    
    template <typename RealType, uint32_t NumSpectralSamples>
    void UpsampledContinuousSpectrumTemplate<RealType, NumSpectralSamples>::computeAdjacents(RealType u, RealType v) {
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
    
    template <typename RealType, uint32_t NumSpectralSamples>
    UpsampledContinuousSpectrumTemplate<RealType, NumSpectralSamples>::UpsampledContinuousSpectrumTemplate(SpectrumType spType, ColorSpace space, RealType e0, RealType e1, RealType e2) {
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
        m_scale = brightness / EqualEnergyReflectance;
        RealType xy[2] = {x, y};
        RealType uv[2];
        xy_to_uv(xy, uv);
        SLRAssert(!std::isinf(uv[0]) && !std::isnan(uv[0]) &&
                  !std::isinf(uv[1]) && !std::isnan(uv[1]) &&
                  !std::isinf(m_scale) && !std::isnan(m_scale), "Invalid value.");
        
        computeAdjacents(uv[0], uv[1]);
    }
    
    template <typename RealType, uint32_t NumSpectralSamples>
    void UpsampledContinuousSpectrumTemplate<RealType, NumSpectralSamples>::calcBounds(uint32_t numBins, RealType* bounds) const {
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
        
        const RealType BinWidth = (WavelengthHighBound - WavelengthLowBound) / numBins;
        const RealType SampleInterval = (MaxWavelength - MinWavelength) / (NumWavelengthSamples - 1);
        
        for (int binIdx = 0; binIdx < numBins; ++binIdx) {
            RealType wlLowInBin = WavelengthLowBound + BinWidth * binIdx;
            RealType wlHighInBin = WavelengthLowBound + BinWidth * (binIdx + 1);
            
            int lowestWlIdxInBin = std::max(0, (int)((wlLowInBin - MinWavelength) / SampleInterval) + 1);
            int highestWlIdxInBin = std::min((int)(NumWavelengthSamples - 1), (int)((wlHighInBin - MinWavelength) / SampleInterval));
            int highestWlIdxOutBin = std::max(0, lowestWlIdxInBin - 1);
            int lowestWlIdxOutBin = std::min((int)(NumWavelengthSamples - 1), highestWlIdxInBin + 1);
            RealType lerpParamAtLowBoundary = (wlLowInBin - (MinWavelength + highestWlIdxOutBin * SampleInterval)) / SampleInterval;
            RealType lerpParamAtHighBoundary = (wlHighInBin - (MinWavelength + highestWlIdxInBin * SampleInterval)) / SampleInterval;
            
            RealType maxValue = -INFINITY;
            RealType lerpedValueAtLow = 0;
            for (int i = 0; i < numAdjacents; ++i) {
                const RealType* spectrum = spectrum_data_points[adjIndices[i]].spectrum;
                lerpedValueAtLow += weights[i] * (spectrum[highestWlIdxOutBin] * (1 - lerpParamAtLowBoundary) + spectrum[lowestWlIdxInBin] * lerpParamAtLowBoundary);
            }
            maxValue = std::max(maxValue, lerpedValueAtLow);
            for (int wlIdx = lowestWlIdxInBin; wlIdx <= highestWlIdxInBin; ++wlIdx) {
                RealType value = 0;
                for (int i = 0; i < numAdjacents; ++i) {
                    const RealType* spectrum = spectrum_data_points[adjIndices[i]].spectrum;
                    lerpedValueAtLow += weights[i] * spectrum[wlIdx];
                }
                maxValue = std::max(maxValue, value);
            }
            RealType lerpedValueAtHigh = 0;
            for (int i = 0; i < numAdjacents; ++i) {
                const RealType* spectrum = spectrum_data_points[adjIndices[i]].spectrum;
                lerpedValueAtHigh += weights[i] * (spectrum[highestWlIdxInBin] * (1 - lerpParamAtHighBoundary) + spectrum[lowestWlIdxOutBin] * lerpParamAtHighBoundary);
            }
            maxValue = std::max(maxValue, lerpedValueAtHigh);
            
            bounds[binIdx] = maxValue;
        }
    }
    
    template <typename RealType, uint32_t NumSpectralSamples>
    SampledSpectrumTemplate<RealType, NumSpectralSamples> 
    UpsampledContinuousSpectrumTemplate<RealType, NumSpectralSamples>::evaluate(const WavelengthSamplesTemplate<RealType, NumSpectralSamples> &wls) const {
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
        
        SampledSpectrumTemplate<RealType, NumSpectralSamples> ret(0.0);
        for (int i = 0; i < WavelengthSamplesTemplate<RealType, NumSpectralSamples>::NumComponents; ++i) {
            RealType lambda = wls[i];
            RealType p = (lambda - MinWavelength) / (MaxWavelength - MinWavelength);
            SLRAssert(p >= 0 && p <= 1, "Wavelength is out of valid range.");
            RealType sBinF = p * (NumWavelengthSamples - 1);
            uint32_t sBin = std::min((uint32_t)sBinF, NumWavelengthSamples - 1);
            uint32_t sBinNext = std::min(sBin + 1, NumWavelengthSamples - 1);
            RealType t = sBinF - sBin;
            for (int j = 0; j < numAdjacents; ++j) {
                const RealType* spectrum = spectrum_data_points[adjIndices[j]].spectrum;
                ret[i] += weights[j] * (spectrum[sBin] * (1 - t) + spectrum[sBinNext] * t);
            }
        }
        
        return ret * m_scale;
    }
    
    template <typename RealType, uint32_t NumSpectralSamples>
    ContinuousSpectrumTemplate<RealType, NumSpectralSamples>* 
    UpsampledContinuousSpectrumTemplate<RealType, NumSpectralSamples>::createScaledAndOffset(RealType scale, RealType offset) const {
        if (scale > 0 && offset == 0)
            return new UpsampledContinuousSpectrumTemplate(m_adjIndices, m_s, m_t, this->m_scale * scale);
        else
            return new ScaledAndOffsetUpsampledContinuousSpectrumTemplate<RealType, NumSpectralSamples>(*this, scale, offset);
    }
    
    template class SLR_API UpsampledContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template class SLR_API UpsampledContinuousSpectrumTemplate<double, NumSpectralSamples>;
    
    
    
    template <typename RealType, uint32_t NumSpectralSamples>
    SampledSpectrumTemplate<RealType, NumSpectralSamples> 
    ScaledAndOffsetUpsampledContinuousSpectrumTemplate<RealType, NumSpectralSamples>::evaluate(const WavelengthSamplesTemplate<RealType, NumSpectralSamples> &wls) const {
        return m_baseSpectrum.evaluate(wls) * m_scale + SampledSpectrumTemplate<RealType, NumSpectralSamples>(m_offset);
    }
    
    template <typename RealType, uint32_t NumSpectralSamples>
    ContinuousSpectrumTemplate<RealType, NumSpectralSamples>* 
    ScaledAndOffsetUpsampledContinuousSpectrumTemplate<RealType, NumSpectralSamples>::createScaledAndOffset(RealType scale, RealType offset) const {
        // (spectrum * scale0 + offset0) * scale1 + offset1 = spectrum * (scale0 * scale1) + (offset0 * scale1 + offset1)
        return new ScaledAndOffsetUpsampledContinuousSpectrumTemplate(m_baseSpectrum, m_scale * scale, m_offset * scale + offset);
    }
    
    template struct SLR_API ScaledAndOffsetUpsampledContinuousSpectrumTemplate<float, NumSpectralSamples>;
    template struct SLR_API ScaledAndOffsetUpsampledContinuousSpectrumTemplate<double, NumSpectralSamples>;
    

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
    
    template <typename RealType, uint32_t NumSpectralSamples>
    SampledSpectrumTemplate<RealType, NumSpectralSamples> sqrt(const SampledSpectrumTemplate<RealType, NumSpectralSamples> &value) {
        SampledSpectrumTemplate<RealType, NumSpectralSamples> ret;
        for (int i = 0; i < NumSpectralSamples; ++i)
            ret[i] = std::sqrt(value[i]);
        return ret;
    }
    template SLR_API SampledSpectrumTemplate<float, NumSpectralSamples> sqrt(const SampledSpectrumTemplate<float, NumSpectralSamples> &value);
    template SLR_API SampledSpectrumTemplate<double, NumSpectralSamples> sqrt(const SampledSpectrumTemplate<double, NumSpectralSamples> &value);
    
    template <typename RealType, uint32_t NumSpectralSamples>
    SampledSpectrumTemplate<RealType, NumSpectralSamples> exp(const SampledSpectrumTemplate<RealType, NumSpectralSamples> &value) {
        SampledSpectrumTemplate<RealType, NumSpectralSamples> ret;
        for (int i = 0; i < NumSpectralSamples; ++i)
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
    
    template class SLR_API SpectrumStorageTemplate<float, NumStrataForStorage>;
    template class SLR_API SpectrumStorageTemplate<double, NumStrataForStorage>;
}
