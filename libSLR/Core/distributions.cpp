//
//  distributions.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "distributions.h"

#include "../BasicTypes/CompensatedSum.h"
#include "../Helper/bmp_exporter.h"
#include "../RNG/LinearCongruentialRNG.h"

namespace SLR {
    template <typename RealType>
    uint32_t sampleDiscrete(const RealType* importances, uint32_t numImportances, RealType u, 
                            RealType* prob, RealType* sumImportances, RealType* remapped) {
        CompensatedSum<RealType> sum(0);
        for (int i = 0; i < numImportances; ++i)
            sum += importances[i];
        *sumImportances = sum;
        
        RealType base = 0;
        RealType su = u * sum;
        CompensatedSum<RealType> cum(0);
        for (int i = 0; i < numImportances; ++i) {
            base = cum;
            cum += importances[i];
            if (su < cum.result) {
                *prob = importances[i] / sum;
                *remapped = (su - base) / importances[i];
                return i;
            }
        }
        *prob = importances[0] / sum;
        return 0;
    }
    template SLR_API uint32_t sampleDiscrete(const float* importances, uint32_t numImportances, float u, 
                                             float* prob, float* sumImportances, float* remapped);
    template SLR_API uint32_t sampleDiscrete(const double* importances, uint32_t numImportances, double u, 
                                             double* prob, double* sumImportances, double* remapped);
    
    
    
    template <typename RealType>
    SLR_API RealType evaluateProbability(const RealType* importances, uint32_t numImportances, uint32_t idx) {
        CompensatedSum<RealType> sum(0);
        for (int i = 0; i < numImportances; ++i)
            sum += importances[i];
        return importances[idx] / sum;
    }
    template SLR_API float evaluateProbability(const float* importances, uint32_t numImportances, uint32_t idx);
    template SLR_API double evaluateProbability(const double* importances, uint32_t numImportances, uint32_t idx);
    
    
    
    // "A Low Distortion Map Between Disk and Square"
    template <typename RealType>
    void concentricSampleDisk(RealType u0, RealType u1, RealType* dx, RealType* dy) {
        RealType r, theta;
        RealType sx = 2 * u0 - 1;
        RealType sy = 2 * u1 - 1;
        
        if (sx == 0 && sy == 0) {
            *dx = 0;
            *dy = 0;
            return;
        }
        if (sx >= -sy) { // region 1 or 2
            if (sx > sy) { // region 1
                r = sx;
                theta = sy / sx;
            }
            else { // region 2
                r = sy;
                theta = 2 - sx / sy;
            }
        }
        else { // region 3 or 4
            if (sx > sy) {/// region 4
                r = -sy;
                theta = 6 + sx / sy;
            }
            else {// region 3
                r = -sx;
                theta = 4 + sy / sx;
            }
        }
        theta *= M_PI_4;
        *dx = r * cos(theta);
        *dy = r * sin(theta);
    }
    template SLR_API void concentricSampleDisk(float u0, float u1, float* dx, float* dy);
    template SLR_API void concentricSampleDisk(double u0, double u1, double* dx, double* dy);
    
    
    
    template <typename RealType>
    DiscreteDistribution1DTemplate<RealType>::DiscreteDistribution1DTemplate(const RealType* values, size_t numValues) {
        m_numValues = (uint32_t)numValues;
        m_PMF = new RealType[m_numValues];
        m_CDF = new RealType[m_numValues + 1];
        std::memcpy(m_PMF, values, sizeof(RealType) * m_numValues);
        
        CompensatedSum<RealType> sum(0);
        m_CDF[0] = 0;
        for (int i = 0; i < m_numValues; ++i) {
            sum += m_PMF[i];
            m_CDF[i + 1] = sum;
        }
        m_integral = sum;
        for (int i = 0; i < m_numValues; ++i) {
            m_PMF[i] /= m_integral;
            m_CDF[i + 1] /= m_integral;
        }
    }
    
    template <typename RealType>
    DiscreteDistribution1DTemplate<RealType>::DiscreteDistribution1DTemplate(const std::vector<RealType> &values) {
        m_numValues = (uint32_t)values.size();
        m_PMF = new RealType[m_numValues];
        m_CDF = new RealType[m_numValues + 1];
        std::memcpy(m_PMF, values.data(), sizeof(RealType) * m_numValues);
        
        CompensatedSum<RealType> sum(0);
        m_CDF[0] = 0;
        for (int i = 0; i < m_numValues; ++i) {
            sum += m_PMF[i];
            m_CDF[i + 1] = sum;
        }
        m_integral = sum;
        for (int i = 0; i < m_numValues; ++i) {
            m_PMF[i] /= m_integral;
            m_CDF[i + 1] /= m_integral;
        }
    }
    
    template <typename RealType>
    uint32_t DiscreteDistribution1DTemplate<RealType>::sample(RealType u, RealType* prob) const {
        SLRAssert(u >= 0 && u < 1, "\"u\" must be in range [0, 1).");
        int idx = m_numValues;
        for (int d = prevPowerOf2(m_numValues); d > 0; d >>= 1)
            if (idx - d > 0 && m_CDF[idx - d] >= u)
                idx -= d;
        --idx;
        *prob = m_PMF[idx];
        return idx;
    };
    
    template <typename RealType>
    uint32_t DiscreteDistribution1DTemplate<RealType>::sample(RealType u, RealType* prob, RealType* remapped) const {
        SLRAssert(u >= 0 && u < 1, "\"u\" must be in range [0, 1).");
        int idx = m_numValues;
        for (int d = prevPowerOf2(m_numValues); d > 0; d >>= 1)
            if (idx - d > 0 && m_CDF[idx - d] >= u)
                idx -= d;
        --idx;
        *prob = m_PMF[idx];
        *remapped = (u - m_CDF[idx]) / (m_CDF[idx + 1] - m_CDF[idx]);
        return idx;
    };
    
    template class SLR_API DiscreteDistribution1DTemplate<float>;
    template class SLR_API DiscreteDistribution1DTemplate<double>;
    
    
    
    template class SLR_API ContinuousDistribution1DTemplate<float>;
    template class SLR_API ContinuousDistribution1DTemplate<double>;
    
    template <typename RealType>
    RegularConstantContinuousDistribution1DTemplate<RealType>::RegularConstantContinuousDistribution1DTemplate(uint32_t numValues, const std::function<RealType(uint32_t)> &pickFunc) :
    m_numValues(numValues) {
        m_PDF = new RealType[m_numValues];
        m_CDF = new RealType[m_numValues + 1];
        for (int i = 0; i < numValues; ++i)
            m_PDF[i] = pickFunc(i);
        
        CompensatedSum<RealType> sum(0);
        m_CDF[0] = 0;
        for (int i = 0; i < m_numValues; ++i) {
            sum += m_PDF[i] / m_numValues;
            m_CDF[i + 1] = sum;
        }
        m_integral = sum;
        if (m_integral > 0) {
            for (int i = 0; i < m_numValues; ++i) {
                m_PDF[i] /= sum;
                m_CDF[i + 1] /= sum;
            }   
        }
    };
    
    template <typename RealType>
    RegularConstantContinuousDistribution1DTemplate<RealType>::RegularConstantContinuousDistribution1DTemplate(const std::vector<RealType> &values) :
    m_numValues((uint32_t)values.size()) {
        m_PDF = new RealType[m_numValues];
        m_CDF = new RealType[m_numValues + 1];
        std::memcpy(m_PDF, values.data(), sizeof(RealType) * m_numValues);
        
        CompensatedSum<RealType> sum{0};
        m_CDF[0] = 0;
        for (int i = 0; i < m_numValues; ++i) {
            sum += m_PDF[i] / m_numValues;
            m_CDF[i + 1] = sum;
        }
        m_integral = sum;
        for (int i = 0; i < m_numValues; ++i) {
            m_PDF[i] /= sum;
            m_CDF[i + 1] /= sum;
        }
    };
    
    template <typename RealType>
    RealType RegularConstantContinuousDistribution1DTemplate<RealType>::sample(RealType u, RealType* PDF) const {
        SLRAssert(u < 1, "\"u\" must be in range [0, 1).");
        int idx = m_numValues;
        for (int d = prevPowerOf2(m_numValues); d > 0; d >>= 1)
            if (idx - d > 0 && m_CDF[idx - d] >= u)
                idx -= d;
        --idx;
        *PDF = m_PDF[idx];
        RealType t = (u - m_CDF[idx]) / (m_CDF[idx + 1] - m_CDF[idx]);
        return (idx + t) / m_numValues;
    };
    
    template <typename RealType>
    RealType RegularConstantContinuousDistribution1DTemplate<RealType>::evaluatePDF(RealType smp) const {
        SLRAssert(smp >= 0 && smp < 1.0, "\"smp\" is out of range [0, 1)");
        return m_PDF[(int32_t)(smp * m_numValues)];
    };
    
    template class SLR_API RegularConstantContinuousDistribution1DTemplate<float>;
    template class SLR_API RegularConstantContinuousDistribution1DTemplate<double>;
    
    
    
    template class SLR_API ContinuousDistribution2DTemplate<float>;
    template class SLR_API ContinuousDistribution2DTemplate<double>;
    
    template <typename RealType>
    RegularConstantContinuousDistribution2DTemplate<RealType>::RegularConstantContinuousDistribution2DTemplate(uint32_t numD1, uint32_t numD2, const std::function<RealType(uint32_t, uint32_t)> &pickFunc) :
    m_num1DDists(numD2) {
        // JP: まず各行に関するDistribution1Dを作成する。
        // EN: First, create Distribution1D's for every rows.
        m_1DDists = (RegularConstantContinuousDistribution1DTemplate<RealType>*)malloc(sizeof(RegularConstantContinuousDistribution1DTemplate<RealType>) * numD2);
        CompensatedSum<RealType> sum(0);
        for (int i = 0; i < numD2; ++i) {
            auto pickFunc1D = std::bind(pickFunc, std::placeholders::_1, i);
            new (m_1DDists + i) RegularConstantContinuousDistribution1DTemplate<RealType>(numD1, pickFunc1D);
            sum += m_1DDists[i].integral();
        }
        m_integral = sum;
        
        // JP: 各行の積分値を用いてDistribution1Dを作成する。
        // EN: create a Distribution1D using integral values of each row.
        auto pickFuncTop = [this](uint32_t idx) { return m_1DDists[idx].integral(); };
        m_top1DDist = new RegularConstantContinuousDistribution1DTemplate<RealType>(numD2, pickFuncTop);
        SLRAssert(std::isfinite(m_integral), "invalid integral value.");
    };
    
    template <typename RealType>
    void RegularConstantContinuousDistribution2DTemplate<RealType>::sample(RealType u0, RealType u1, RealType* d0, RealType* d1, RealType* PDF) const {
        SLRAssert(u0 >= 0 && u0 < 1, "\"u0\" must be in range [0, 1).");
        SLRAssert(u1 >= 0 && u1 < 1, "\"u1\" must be in range [0, 1).");
        RealType topPDF;
        *d1 = m_top1DDist->sample(u1, &topPDF);
        uint32_t idx1D = std::min(uint32_t(m_num1DDists * *d1), m_num1DDists - 1);
        *d0 = m_1DDists[idx1D].sample(u0, PDF);
        *PDF *= topPDF;
    };
    
    template <typename RealType>
    RealType RegularConstantContinuousDistribution2DTemplate<RealType>::evaluatePDF(RealType d0, RealType d1) const {
        SLRAssert(d0 >= 0 && d0 < 1.0, "\"d0\" is out of range [0, 1)");
        SLRAssert(d1 >= 0 && d1 < 1.0, "\"d1\" is out of range [0, 1)");
        uint32_t idx1D = std::min(uint32_t(m_num1DDists * d1), m_num1DDists - 1);
        return m_top1DDist->evaluatePDF(d1) * m_1DDists[idx1D].evaluatePDF(d0);
    };
    
    // For debug visualization.
    template <typename RealType>
    void RegularConstantContinuousDistribution2DTemplate<RealType>::exportBMP(const std::string &filename, bool logScale, float gamma) const {
        uint32_t width = m_1DDists[0].numValues();
        uint32_t height = m_num1DDists;
        uint32_t byteWidth = width * 3 + width % 4;
        uint8_t* data = (uint8_t*)malloc(height * byteWidth);
        
        float minValue = INFINITY;
        float maxValue = -INFINITY;
        for (int i = 0; i < height; ++i) {
            const RealType* PDF = m_1DDists[i].PDF();
            float ratio = m_top1DDist->PDF()[i];
            for (int j = 0; j < width; ++j) {
                float value = ratio * PDF[j];
                if (logScale)
                    value = std::log(value);
                if (std::isfinite(value)) {
                    minValue = std::min(minValue, value);
                    maxValue = std::max(maxValue, value);   
                }
            }
        }
        for (int i = 0; i < height; ++i) {
            const RealType* PDF = m_1DDists[i].PDF();
            float ratio = m_top1DDist->PDF()[i];
            for (int j = 0; j < width; ++j) {
                float value = ratio * PDF[j];
                if (logScale) {
                    value = std::log(value);
                    if (std::isfinite(value)) {
                        value = (value - minValue) / (maxValue - minValue);
                    }
                    else {
                        value = 0;
                    }
                }
                else {
                    value /= maxValue;
                }
                uint8_t pixVal = uint8_t(std::pow(value, 1.0f / gamma) * 255);
                
                uint32_t idx = (height - i - 1) * byteWidth + 3 * j;
                data[idx + 0] = pixVal;
                data[idx + 1] = pixVal;
                data[idx + 2] = pixVal;
            }
        }
        saveBMP(filename.c_str(), data, width, height);
        free(data);
    };
    
    template class SLR_API RegularConstantContinuousDistribution2DTemplate<float>;
    template class SLR_API RegularConstantContinuousDistribution2DTemplate<double>;
    
    template <typename RealType>
    MultiContinuousDistribution2DTemplate<RealType>::MultiContinuousDistribution2DTemplate(const ContinuousDistribution2DTemplate<RealType>** dists, const RealType* importances, uint32_t numDists) : 
    m_selectDist(DiscreteDistribution1DTemplate<RealType>(importances, numDists)) {
        SLRAssert(numDists <= MaxNumDist2Ds, "Number of distributions is limited to %u", MaxNumDist2Ds);
        for (int i = 0; i < numDists; ++i)
            m_dist2Ds[i] = dists[i];
        for (int i = numDists; i < MaxNumDist2Ds; ++i)
            m_dist2Ds[i] = nullptr;
    }
    
    template <typename RealType>
    void MultiContinuousDistribution2DTemplate<RealType>::sample(RealType u0, RealType u1, RealType* d0, RealType* d1, RealType* PDF) const {
        SLRAssert(u0 >= 0 && u0 < 1, "\"u0\" must be in range [0, 1).");
        SLRAssert(u1 >= 0 && u1 < 1, "\"u1\" must be in range [0, 1).");
        RealType selectProb;
        uint32_t selIdx = m_selectDist.sample(u0, &selectProb, &u0); 
        
        m_dist2Ds[selIdx]->sample(u0, u1, d0, d1, PDF);
        *PDF *= selectProb;
        
        for (int i = 0; i < m_selectDist.numValues(); ++i) {
            if (i == selIdx)
                continue;
            *PDF += m_selectDist.evaluatePMF(i) * m_dist2Ds[i]->evaluatePDF(*d0, *d1);
        }
    };
    
    template <typename RealType>
    RealType MultiContinuousDistribution2DTemplate<RealType>::evaluatePDF(RealType d0, RealType d1) const {
        SLRAssert(d0 >= 0 && d0 < 1.0, "\"d0\" is out of range [0, 1)");
        SLRAssert(d1 >= 0 && d1 < 1.0, "\"d1\" is out of range [0, 1)");
        
        RealType retPDF = 0;
        for (int i = 0; i < m_selectDist.numValues(); ++i)
            retPDF += m_selectDist.evaluatePMF(i) * m_dist2Ds[i]->evaluatePDF(d0, d1);
        
        return retPDF;
    };
    
    template class SLR_API MultiContinuousDistribution2DTemplate<float>;
    template class SLR_API MultiContinuousDistribution2DTemplate<double>;
    
    
    
    template <typename RealType>
    RealType ImprovedPerlinNoise3DGeneratorTemplate<RealType>::gradient(uint32_t hash, RealType xu, RealType yu, RealType zu) {
        switch (hash & 0xF) {
            case 0x0: return  xu + yu;
            case 0x1: return -xu + yu;
            case 0x2: return  xu - yu;
            case 0x3: return -xu - yu;
            case 0x4: return  xu + zu;
            case 0x5: return -xu + zu;
            case 0x6: return  xu - zu;
            case 0x7: return -xu - zu;
            case 0x8: return  yu + zu;
            case 0x9: return -yu + zu;
            case 0xA: return  yu - zu;
            case 0xB: return -yu - zu;
            case 0xC: return  yu + xu;
            case 0xD: return -yu + zu;
            case 0xE: return  yu - xu;
            case 0xF: return -yu - zu;
            default: return 0; // never happens
        }
    }
    
    template <typename RealType>
    RealType ImprovedPerlinNoise3DGeneratorTemplate<RealType>::evaluate(const Point3DTemplate<RealType> &p) const {        
        RealType x = p.x, y = p.y, z = p.z;
        
        // If we have any repeat on, change the coordinates to their "local" repetitions.
        if (m_repeat > 0) {
            x = std::fmod(x, m_repeat);
            y = std::fmod(y, m_repeat);
            z = std::fmod(z, m_repeat);
            if (x < 0)
                x += m_repeat;
            if (y < 0)
                y += m_repeat;
            if (z < 0)
                z += m_repeat;
        }
        
        int32_t iEvalCoord[3];
        iEvalCoord[0] = std::floor(x);
        iEvalCoord[1] = std::floor(y);
        iEvalCoord[2] = std::floor(z);
        
        const auto fade = [](RealType t) {
            // Fade function as defined by Ken Perlin.
            // This eases coordinate values so that they will "ease" towards integral values.
            // This ends up smoothing the final output.
            // 6t^5 - 15t^4 + 10t^3
            return t * t * t * (t * (t * 6 - 15) + 10);
        };
        
        // Calculate the "unit cube" that the point asked will be located in.
        // The left bound is ( |_x_|,|_y_|,|_z_| ) and the right bound is that plus 1. 
        // Next we calculate the location (from 0.0 to 1.0) in that cube.
        // We also fade the location to smooth the result.
        RealType xu = x - iEvalCoord[0];
        RealType yu = y - iEvalCoord[1];
        RealType zu = z - iEvalCoord[2];
        SLRAssert(xu >= 0 && xu <= 1 && yu >= 0 && yu <= 1 && zu >= 0 && zu <= 1, "xu, yu, zu must be in the unit cube [0, 1]^3.");
        RealType u = fade(xu);
        RealType v = fade(yu);
        RealType w = fade(zu);
        
        // Code in the reference site uses a permutation table of size 256.
        // It imposes the repetition pattern of size 256 even if repetition is not specified.
        // Modified it by instead using a pseudo random number generator based on integer coordinate.
        uint8_t cornerHashes[8];
        for (int iz = 0; iz < 2; ++iz) {
            for (int iy = 0; iy < 2; ++iy) {
                for (int ix = 0; ix < 2; ++ix) {
                    int32_t iCoord[3] = {iEvalCoord[0] + ix, iEvalCoord[1] + iy, iEvalCoord[2] + iz};
                    if (m_repeat > 0) {
                        iCoord[0] %= m_repeat;
                        iCoord[1] %= m_repeat;
                        iCoord[2] %= m_repeat;
                    }
                    uint32_t hash = getFNV1Hash32((uint8_t*)iCoord, sizeof(iCoord));
                    LinearCongruentialRNG rng(hash);
                    cornerHashes[4 * iz + 2 * iy + ix] = (rng.getUInt() >> 24) & 0xF;
                }
            }
        }
        
        const auto lerp = [](RealType v0, RealType v1, RealType t) {
            return v0 * (1 - t) + v1 * t;
        };
        
        // The gradient function calculates the dot product between a pseudorandom gradient vector and 
        // the vector from the input coordinate to the 8 surrounding points in its unit cube.
        // This is all then lerped together as a sort of weighted average based on the faded (u,v,w) values we made earlier.
        RealType _llValue = lerp(gradient(cornerHashes[0], xu, yu, zu), gradient(cornerHashes[1], xu - 1, yu, zu), u);
        RealType _ulValue = lerp(gradient(cornerHashes[2], xu, yu - 1, zu), gradient(cornerHashes[3], xu - 1, yu - 1, zu), u);
        RealType __lValue = lerp(_llValue, _ulValue, v);
        
        RealType _luValue = lerp(gradient(cornerHashes[4], xu, yu, zu - 1), gradient(cornerHashes[5], xu - 1, yu, zu - 1), u);
        RealType _uuValue = lerp(gradient(cornerHashes[6], xu, yu - 1, zu - 1), gradient(cornerHashes[7], xu - 1, yu - 1, zu - 1), u);
        RealType __uValue = lerp(_luValue, _uuValue, v);
        
        // For convenience we bound it to 0 - 1 (theoretical min/max before is -1 - 1)
        RealType ret = (lerp(__lValue, __uValue, w) + 1) / 2;
        SLRAssert(ret >= 0 && ret <= 1.0f, "Return value is invalid.");
        return ret;
    }
    
    template class SLR_API ImprovedPerlinNoise3DGeneratorTemplate<float>;
    template class SLR_API ImprovedPerlinNoise3DGeneratorTemplate<double>;
    
    
    
    template <typename RealType>
    RealType MultiOctavePerlinNoise3DGeneratorTemplate<RealType>::evaluate(const Point3DTemplate<RealType> &p) const {
        RealType total = 0;
        RealType frequency = m_initialFrequency;
        RealType amplitude = m_initialAmplitude;
        for (int i = 0; i < m_numOctaves; ++i) {
            total += m_primaryNoiseGen.evaluate(frequency * p) * amplitude;
            
            amplitude *= m_persistence;
            frequency *= m_frequencyMultiplier;
        }
        
        return total;
    }
    
    template class SLR_API MultiOctavePerlinNoise3DGeneratorTemplate<float>;
    template class SLR_API MultiOctavePerlinNoise3DGeneratorTemplate<double>;
    
    
    
    template <typename RealType>
    Vector3DTemplate<RealType> CurlNoise3DGeneratorTemplate<RealType>::evaluate(const Point3DTemplate<RealType> &p) const {
        using Point3DType = Point3DTemplate<RealType>;
        using Vector3DType = Vector3DTemplate<RealType>;
        
        const RealType Offset = 50;
        const RealType Delta = 1e-3 / m_maxFrequency;
        
        const RealType Psi1Base = m_noiseGen.evaluate(Point3DType(p.x, p.y, p.z));
        const RealType Psi2Base = m_noiseGen.evaluate(Point3DType(p.x + Offset, p.y, p.z));
        const RealType Psi3Base = m_noiseGen.evaluate(Point3DType(p.x, p.y + Offset, p.z));
        
        RealType rPsi1ry = (m_noiseGen.evaluate(Point3DType(p.x, p.y + Delta, p.z)) - Psi1Base) / Delta;
        RealType rPsi1rz = (m_noiseGen.evaluate(Point3DType(p.x, p.y, p.z + Delta)) - Psi1Base) / Delta;
        RealType rPsi2rx = (m_noiseGen.evaluate(Point3DType(p.x + Offset + Delta, p.y, p.z)) - Psi2Base) / Delta;
        RealType rPsi2rz = (m_noiseGen.evaluate(Point3DType(p.x + Offset, p.y, p.z + Delta)) - Psi2Base) / Delta;
        RealType rPsi3rx = (m_noiseGen.evaluate(Point3DType(p.x + Offset + Delta, p.y, p.z)) - Psi3Base) / Delta;
        RealType rPsi3ry = (m_noiseGen.evaluate(Point3DType(p.x + Offset, p.y + Delta, p.z)) - Psi3Base) / Delta;
        
        return Vector3DType(rPsi3ry - rPsi2rz, rPsi1rz - rPsi3rx, rPsi2rx - rPsi1ry);
    }
    
    template class SLR_API CurlNoise3DGeneratorTemplate<float>;
    template class SLR_API CurlNoise3DGeneratorTemplate<double>;
    
    
    
    template <typename RealType>
    void WorleyNoise3DGeneratorTemplate<RealType>::evaluate(const Point3DTemplate<RealType> &p, RealType* closestDistance, uint32_t* hashOfClosest, uint32_t* closestFPIdx) const {
        using Point3DType = Point3DTemplate<RealType>;
        
        int32_t iEvalCoord[3];
        iEvalCoord[0] = std::floor(p.x);
        iEvalCoord[1] = std::floor(p.y);
        iEvalCoord[2] = std::floor(p.z);
        
        // JP: 本来は5 x 5 x 5 - 8 = 127の隣接セルを評価する必要がある。
        //     評価点があるセル内において評価点が8区画のどれに属するかの情報を用いて最適化しても56セルは評価する必要がある。
        //     が、27セルでもほとんどうまくいく。
        // EN: Strictly, this requires evaluation of 5 x 5 x 5 - 8 = 127 adjacent cells.
        //     It is necessary to evaluate 56 cells even if 
        //     optimization is applied using information that which of 8 blocks the evaluation point belongs to in the cell where the evaluation point is in.  
        //     However, it seems evaluating only 27 cells works well.
        *closestDistance = INFINITY;
        for (int iz = -1; iz <= 1; ++iz) {
            for (int iy = -1; iy <= 1; ++iy) {
                for (int ix = -1; ix <= 1; ++ix) {
                    int32_t iCoord[3] = {iEvalCoord[0] + ix, iEvalCoord[1] + iy, iEvalCoord[2] + iz};
                    uint32_t hash = getFNV1Hash32((uint8_t*)iCoord, sizeof(iCoord));
                    LinearCongruentialRNG rng(hash);
                    
                    uint32_t numFeaturePoints = 1 + std::min(int32_t(8 * rng.getFloat0cTo1o()), 8);
                    for (int i = 0; i < numFeaturePoints; ++i) {
                        Point3DType fp = Point3DType(iCoord[0] + rng.getFloat0cTo1o(), 
                                                     iCoord[1] + rng.getFloat0cTo1o(), 
                                                     iCoord[2] + rng.getFloat0cTo1o());
                        RealType dist = distance(p, fp);
                        if (dist < *closestDistance) {
                            *closestDistance = dist;
                            *hashOfClosest = hash;
                            *closestFPIdx = i;
                        }
                    }
                }
            }
        }
    }
    
    template class SLR_API WorleyNoise3DGeneratorTemplate<float>;
    template class SLR_API WorleyNoise3DGeneratorTemplate<double>;
    
    
    
    template <typename RealType>
    RealType MultiOctaveWorleyNoise3DGeneratorTemplate<RealType>::evaluate(const Point3DTemplate<RealType> &p) const {
        RealType total = 0;
        RealType frequency = m_initialFrequency;
        RealType amplitude = m_initialAmplitude;
        for (int i = 0; i < m_numOctaves; ++i) {
            RealType closestDistance;
            uint32_t hashOfClosest;
            uint32_t closestFPIdx;
            m_primaryNoiseGen.evaluate(frequency * p, &closestDistance, &hashOfClosest, &closestFPIdx);
            total += (closestDistance / std::sqrt(3.0)) * amplitude;
            
            amplitude *= m_persistence;
            frequency *= m_frequencyMultiplier;
        }
        
        return std::fmin(total, m_clipValue);
    }
    
    template class SLR_API MultiOctaveWorleyNoise3DGeneratorTemplate<float>;
    template class SLR_API MultiOctaveWorleyNoise3DGeneratorTemplate<double>;
}
