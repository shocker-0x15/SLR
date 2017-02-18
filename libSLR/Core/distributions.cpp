//
//  distributions.cpp
//
//  Created by 渡部 心 on 2016/01/26.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "distributions.h"
#include "../BasicTypes/CompensatedSum.h"
#include "../Helper/bmp_exporter.h"

namespace SLR {
    template <typename RealType>
    uint32_t sampleDiscrete(const RealType* importances, RealType* sumImportances, RealType* base, uint32_t numImportances, RealType u) {
        CompensatedSum<RealType> sum(0);
        for (int i = 0; i < numImportances; ++i)
            sum += importances[i];
        *sumImportances = sum.result;
        
        RealType su = u * sum.result;
        CompensatedSum<RealType> cum(0);
        for (int i = 0; i < numImportances; ++i) {
            *base = cum.result;
            cum += importances[i];
            if (su < cum.result)
                return i;
        }
        return 0;
    }
    template SLR_API uint32_t sampleDiscrete(const float* importances, float* sumImportances, float* base, uint32_t numImportances, float u);
    template SLR_API uint32_t sampleDiscrete(const double* importances, double* sumImportances, double* base, uint32_t numImportances, double u);
    
    
    
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
    };
    
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
        for (int i = 0; i < m_numValues; ++i) {
            m_PDF[i] /= sum;
            m_CDF[i + 1] /= sum;
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
    
    
    
    template <typename RealType>
    RegularConstantContinuousDistribution2DTemplate<RealType>::RegularConstantContinuousDistribution2DTemplate(uint32_t numD1, uint32_t numD2, const std::function<RealType(uint32_t, uint32_t)> &pickFunc) :
    m_num1DDists(numD2) {
        m_1DDists = (RegularConstantContinuousDistribution1DTemplate<RealType>*)malloc(sizeof(RegularConstantContinuousDistribution1DTemplate<RealType>) * numD2);
        CompensatedSum<RealType> sum(0);
        for (int i = 0; i < numD2; ++i) {
            auto pickFunc1D = std::bind(pickFunc, std::placeholders::_1, i);
            new (m_1DDists + i) RegularConstantContinuousDistribution1DTemplate<RealType>(numD1, pickFunc1D);
            sum += m_1DDists[i].integral();
        }
        m_integral = sum;
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
    
    template <typename RealType>
    void RegularConstantContinuousDistribution2DTemplate<RealType>::exportBMP(const std::string &filename, float gamma) const {
        uint32_t width = m_1DDists[0].numValues();
        uint32_t height = m_num1DDists;
        uint32_t byteWidth = width * 3 + width % 4;
        uint8_t* data = (uint8_t*)malloc(height * byteWidth);
        
        float maxValue = -INFINITY;
        for (int i = 0; i < height; ++i) {
            const RealType* PDF = m_1DDists[i].PDF();
            float ratio = m_top1DDist->PDF()[i];
            for (int j = 0; j < width; ++j) {
                float value = ratio * PDF[j];
                if (value > maxValue)
                    maxValue = value;
            }
        }
        for (int i = 0; i < height; ++i) {
            const RealType* PDF = m_1DDists[i].PDF();
            float ratio = m_top1DDist->PDF()[i];
            for (int j = 0; j < width; ++j) {
                float value = std::pow(ratio * PDF[j] / maxValue, 1.0f / gamma);
                uint8_t pixVal = uint8_t(value * 255);
                
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
}
