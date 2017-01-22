//
//  distributions.h
//
//  Created by 渡部 心 on 2015/06/28.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__distributions__
#define __SLR__distributions__

#include "../defines.h"
#include "../references.h"

namespace SLR {
    template <typename RealType>
    SLR_API uint32_t sampleDiscrete(const RealType* importances, RealType* sumImportances, RealType* base, uint32_t numImportances, RealType u);
    
    template <typename RealType>
    SLR_API RealType evaluateProbability(const RealType* importances, uint32_t numImportances, uint32_t idx);
    
    template <typename RealType>
    SLR_API void concentricSampleDisk(RealType u0, RealType u1, RealType* dx, RealType* dy);
    
    template <typename RealType>
    inline Vector3Template<RealType> cosineSampleHemisphere(RealType u0, RealType u1) {
//        RealType phi = 2 * M_PI * u1;
//        RealType theta = std::asin(std::sqrt(u0));
//        return Vector3Template<RealType>(std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta));
        RealType x, y;
        concentricSampleDisk(u0, u1, &x, &y);
        return Vector3Template<RealType>(x, y, std::sqrt(std::fmax(0.0f, 1.0f - x * x - y * y)));
    }
    
    template <typename RealType, int N>
    inline Vector3Template<RealType> cosNSampleHemisphere(RealType u0, RealType u1) {
        RealType phi = 2 * M_PI * u1;
        RealType theta = std::acos(std::pow(u0, 1.0 / (1 + N)));
        return Vector3Template<RealType>(std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta));
    }
    
    template <typename RealType>
    inline Vector3Template<RealType> uniformSampleHemisphere(RealType u0, RealType u1) {
        RealType phi = 2 * M_PI * u1;
        RealType theta = std::acos(1 - u0);
        return Vector3Template<RealType>(std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta));
    }
    
    template <typename RealType>
    inline Vector3Template<RealType> uniformSampleSphere(RealType u0, RealType u1) {
        RealType phi = 2 * M_PI * u1;
        RealType theta = std::acos(1 - 2 * u0);
        return Vector3Template<RealType>(std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta));
    }
    
    template <typename RealType>
    inline Vector3Template<RealType> uniformSampleCone(RealType u0, RealType u1, RealType cosThetaMax) {
        RealType phi = 2 * M_PI * u1;
        RealType theta = std::acos(1 - (1 - cosThetaMax) * u0);
        return Vector3Template<RealType>(std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta));
    }
    
    template <typename RealType>
    inline void uniformSampleTriangle(RealType u0, RealType u1, RealType* b0, RealType* b1) {
        RealType su1 = std::sqrt(u0);
        *b0 = 1.0f - su1;
        *b1 = u1 * su1;
    }
    
    
    
    template <typename RealType>
    class SLR_API RegularConstantDiscrete1DTemplate {
        RealType* m_PMF;
        RealType* m_CDF;
        RealType m_integral;
        uint32_t m_numValues;
    public:
        RegularConstantDiscrete1DTemplate(const std::vector<RealType> &values);
        ~RegularConstantDiscrete1DTemplate() {
            delete[] m_PMF;
            delete[] m_CDF;
        };
        
        uint32_t sample(RealType u, RealType* prob) const;
        uint32_t sample(RealType u, RealType* prob, RealType* remapped) const;
        RealType evaluatePMF(uint32_t idx) const {
            SLRAssert(idx >= 0 && idx < m_numValues, "\"idx\" is out of range [0, %u)", m_numValues);
            return m_PMF[idx];
        };
        
        RealType integral() const { return m_integral; };
    };
    
    
    
    template <typename RealType>
    class SLR_API RegularConstantContinuous1DTemplate {
        RealType* m_PDF;
        RealType* m_CDF;
        RealType m_integral;
        uint32_t m_numValues;
    public:
        RegularConstantContinuous1DTemplate(uint32_t numValues, const std::function<RealType(uint32_t)> &pickFunc);
        RegularConstantContinuous1DTemplate(const std::vector<RealType> &values);
        ~RegularConstantContinuous1DTemplate() {
            delete[] m_PDF;
            delete[] m_CDF;
        };
        
        RealType sample(RealType u, RealType* PDF) const;
        RealType evaluatePDF(RealType smp) const;
        RealType integral() const { return m_integral; };
        uint32_t numValues() const { return m_numValues; };
        const RealType* PDF() const { return m_PDF; };
    };
    
    
    
    template <typename RealType>
    class SLR_API RegularConstantContinuous2DTemplate {
        RegularConstantContinuous1DTemplate<RealType>* m_1DDists;
        uint32_t m_num1DDists;
        RealType m_integral;
        RegularConstantContinuous1DTemplate<RealType>* m_top1DDist;
    public:
        RegularConstantContinuous2DTemplate(uint32_t numD1, uint32_t numD2, const std::function<RealType(uint32_t, uint32_t)> &pickFunc);;
        ~RegularConstantContinuous2DTemplate() {
            free(m_1DDists);
            delete m_top1DDist;
        };
        
        void sample(RealType u0, RealType u1, RealType* d0, RealType* d1, RealType* PDF) const;
        RealType evaluatePDF(RealType d0, RealType d1) const;
        RealType integral() const { return m_integral; };
        void exportBMP(const std::string &filename, float gamma = 1.0f) const;
    };    
}

#endif
