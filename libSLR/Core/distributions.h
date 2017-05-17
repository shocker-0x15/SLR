//
//  distributions.h
//
//  Created by 渡部 心 on 2015/06/28.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_distributions__
#define __SLR_distributions__

#include "../defines.h"
#include "../declarations.h"

namespace SLR {
    template <typename RealType>
    SLR_API uint32_t sampleDiscrete(const RealType* importances, uint32_t numImportances, RealType u, 
                                    RealType* prob, RealType* sumImportances, RealType* remapped);
    
    template <typename RealType>
    SLR_API RealType evaluateProbability(const RealType* importances, uint32_t numImportances, uint32_t idx);
    
    
    
    template <typename RealType>
    SLR_API void concentricSampleDisk(RealType u0, RealType u1, RealType* dx, RealType* dy);
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> cosineSampleHemisphere(RealType u0, RealType u1) {
//        RealType phi = 2 * M_PI * u1;
//        RealType theta = std::asin(std::sqrt(u0));
//        return Vector3DTemplate<RealType>(std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta));
        RealType x, y;
        concentricSampleDisk(u0, u1, &x, &y);
        return Vector3DTemplate<RealType>(x, y, std::sqrt(std::fmax(0.0f, 1.0f - x * x - y * y)));
    }
    
    template <typename RealType, int N>
    inline Vector3DTemplate<RealType> cosNSampleHemisphere(RealType u0, RealType u1) {
        RealType phi = 2 * M_PI * u1;
        RealType theta = std::acos(std::pow(u0, 1.0 / (1 + N)));
        return Vector3DTemplate<RealType>(std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta));
    }
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> uniformSampleHemisphere(RealType u0, RealType u1) {
        RealType phi = 2 * M_PI * u1;
        RealType theta = std::acos(1 - u0);
        return Vector3DTemplate<RealType>(std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta));
    }
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> uniformSampleSphere(RealType u0, RealType u1) {
        RealType phi = 2 * M_PI * u1;
        RealType theta = std::acos(1 - 2 * u0);
        return Vector3DTemplate<RealType>(std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta));
    }
    
    template <typename RealType>
    inline Vector3DTemplate<RealType> uniformSampleCone(RealType u0, RealType u1, RealType cosThetaMax) {
        RealType phi = 2 * M_PI * u1;
        RealType theta = std::acos(1 - (1 - cosThetaMax) * u0);
        return Vector3DTemplate<RealType>(std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta));
    }
    
    template <typename RealType>
    inline void uniformSampleTriangle(RealType u0, RealType u1, RealType* b0, RealType* b1) {
        RealType su1 = std::sqrt(u0);
        *b0 = 1.0f - su1;
        *b1 = u1 * su1;
    }
    
    
    
    template <typename RealType>
    class SLR_API DiscreteDistribution1DTemplate {
        RealType* m_PMF;
        RealType* m_CDF;
        RealType m_integral;
        uint32_t m_numValues;
    public:
        DiscreteDistribution1DTemplate() : m_PMF(nullptr), m_CDF(nullptr) { }
        DiscreteDistribution1DTemplate(const RealType* values, size_t numValues);
        DiscreteDistribution1DTemplate(const std::vector<RealType> &values);
        ~DiscreteDistribution1DTemplate() {
            if (m_PMF)
                delete[] m_PMF;
            if (m_CDF)
                delete[] m_CDF;
        }
        
        uint32_t sample(RealType u, RealType* prob) const;
        uint32_t sample(RealType u, RealType* prob, RealType* remapped) const;
        RealType evaluatePMF(uint32_t idx) const {
            SLRAssert(idx >= 0 && idx < m_numValues, "\"idx\" is out of range [0, %u)", m_numValues);
            return m_PMF[idx];
        }
        
        RealType integral() const { return m_integral; }
        uint32_t numValues() const { return m_numValues; }
    };
    
    
    
    template <typename RealType>
    class SLR_API ContinuousDistribution1DTemplate {
    public:
        virtual ~ContinuousDistribution1DTemplate() { }
        virtual RealType sample(RealType u, RealType* PDF) const = 0;
        virtual RealType evaluatePDF(RealType smp) const = 0;
        virtual RealType integral() const = 0;
    };
    
    template <typename RealType>
    class SLR_API RegularConstantContinuousDistribution1DTemplate : public ContinuousDistribution1DTemplate<RealType> {
        RealType* m_PDF;
        RealType* m_CDF;
        RealType m_integral;
        uint32_t m_numValues;
    public:
        RegularConstantContinuousDistribution1DTemplate(uint32_t numValues, const std::function<RealType(uint32_t)> &pickFunc);
        RegularConstantContinuousDistribution1DTemplate(const std::vector<RealType> &values);
        ~RegularConstantContinuousDistribution1DTemplate() {
            delete[] m_PDF;
            delete[] m_CDF;
        }
        
        RealType sample(RealType u, RealType* PDF) const override;
        RealType evaluatePDF(RealType smp) const override;
        RealType integral() const override { return m_integral; }
        
        uint32_t numValues() const { return m_numValues; }
        const RealType* PDF() const { return m_PDF; }
    };
    
    
    
    template <typename RealType>
    class SLR_API ContinuousDistribution2DTemplate {
    public:
        virtual ~ContinuousDistribution2DTemplate() { }
        virtual void sample(RealType u0, RealType u1, RealType* d0, RealType* d1, RealType* PDF) const = 0;
        virtual RealType evaluatePDF(RealType d0, RealType d1) const = 0;
    };
    
    template <typename RealType>
    class SLR_API RegularConstantContinuousDistribution2DTemplate : public ContinuousDistribution2DTemplate<RealType> {
        RegularConstantContinuousDistribution1DTemplate<RealType>* m_1DDists;
        uint32_t m_num1DDists;
        RealType m_integral;
        RegularConstantContinuousDistribution1DTemplate<RealType>* m_top1DDist;
    public:
        RegularConstantContinuousDistribution2DTemplate(uint32_t numD1, uint32_t numD2, const std::function<RealType(uint32_t, uint32_t)> &pickFunc);
        ~RegularConstantContinuousDistribution2DTemplate() {
            free(m_1DDists);
            delete m_top1DDist;
        }
        
        void sample(RealType u0, RealType u1, RealType* d0, RealType* d1, RealType* PDF) const override;
        RealType evaluatePDF(RealType d0, RealType d1) const override;

        void exportBMP(const std::string &filename, bool logScale = false, float gamma = 1.0f) const;
    };
    
    template <typename RealType>
    class SLR_API MultiContinuousDistribution2DTemplate : public ContinuousDistribution2DTemplate<RealType> {
        static const uint32_t MaxNumDist2Ds = 4;
        DiscreteDistribution1DTemplate<RealType> m_selectDist;
        const ContinuousDistribution2DTemplate<RealType>* m_dist2Ds[MaxNumDist2Ds];
    public:
        MultiContinuousDistribution2DTemplate(const ContinuousDistribution2DTemplate<RealType>** dists, const RealType* importances, uint32_t numDists);
        
        void sample(RealType u0, RealType u1, RealType* d0, RealType* d1, RealType* PDF) const override;
        RealType evaluatePDF(RealType d0, RealType d1) const override;
    };
    
    
    
    // This code is based on the web site: adrian's soapbox
    // http://flafla2.github.io/2014/08/09/perlinnoise.html
    template <typename RealType>
    class MultiOctaveImprovedPerlinNoise3DGenerator {
        static const uint8_t PermutationTable[];
        
        static RealType gradient(uint32_t hash, RealType xu, RealType yu, RealType zu);
        
        RealType primaryPerlinNoise(RealType x, RealType y, RealType z) const;
        
        uint32_t m_numOctaves;
        RealType m_initialFrequency;
        RealType m_initialAmplitude;
        RealType m_frequencyMultiplier;
        RealType m_persistence;
        int32_t m_repeat;
        
    public:
        MultiOctaveImprovedPerlinNoise3DGenerator(uint32_t numOctaves, RealType initialFrequency, RealType supValueOrInitialAmplitude, bool supSpecified,  
                                                  RealType frequencyMultiplier, RealType persistence, int32_t repeat) : 
        m_numOctaves(numOctaves), 
        m_initialFrequency(initialFrequency), 
        m_frequencyMultiplier(frequencyMultiplier), m_persistence(persistence), m_repeat(repeat) {
            if (supSpecified) {
                RealType amplitude = 1.0f;
                RealType supValue = 0; // used for normalizing result to 0.0 - 1.0
                for (int i = 0; i < m_numOctaves; ++i) {                
                    supValue += amplitude;
                    amplitude *= m_persistence;
                }
                m_initialAmplitude = supValueOrInitialAmplitude / supValue;   
            }
            else {
                m_initialAmplitude = supValueOrInitialAmplitude;
            }
        }
        
        RealType evaluate(RealType x, RealType y, RealType z) const;
    };
    
    template <typename RealType>
    const uint8_t MultiOctaveImprovedPerlinNoise3DGenerator<RealType>::PermutationTable[] = {
        151, 160, 137,  91,  90,  15, 131,  13, 201,  95,  96,  53, 194, 233,   7, 225,
        140,  36, 103,  30,  69, 142,   8,  99,  37, 240,  21,  10,  23, 190,   6, 148,
        247, 120, 234,  75,   0,  26, 197,  62,  94, 252, 219, 203, 117,  35,  11,  32,
         57, 177,  33,  88, 237, 149,  56,  87, 174,  20, 125, 136, 171, 168,  68, 175,
         74, 165,  71, 134, 139,  48,  27, 166,  77, 146, 158, 231,  83, 111, 229, 122,
         60, 211, 133, 230, 220, 105,  92,  41,  55,  46, 245,  40, 244, 102, 143,  54,
         65,  25,  63, 161,   1, 216,  80,  73, 209,  76, 132, 187, 208,  89,  18, 169,
        200, 196, 135, 130, 116, 188, 159,  86, 164, 100, 109, 198, 173, 186,   3,  64,
         52, 217, 226, 250, 124, 123,   5, 202,  38, 147, 118, 126, 255,  82,  85, 212,
        207, 206,  59, 227,  47,  16,  58,  17, 182, 189,  28,  42, 223, 183, 170, 213,
        119, 248, 152,   2,  44, 154, 163,  70, 221, 153, 101, 155, 167,  43, 172,   9,
        129,  22,  39, 253,  19,  98, 108, 110,  79, 113, 224, 232, 178, 185, 112, 104,
        218, 246,  97, 228, 251,  34, 242, 193, 238, 210, 144,  12, 191, 179, 162, 241,
         81,  51, 145, 235, 249,  14, 239, 107,  49, 192, 214,  31, 181, 199, 106, 157,
        184,  84, 204, 176, 115, 121,  50,  45, 127,   4, 150, 254, 138, 236, 205,  93,
        222, 114,  67,  29,  24,  72, 243, 141, 128, 195,  78,  66, 215,  61, 156, 180,
        
        151, 160, 137,  91,  90,  15, 131,  13, 201,  95,  96,  53, 194, 233,   7, 225,
        140,  36, 103,  30,  69, 142,   8,  99,  37, 240,  21,  10,  23, 190,   6, 148,
        247, 120, 234,  75,   0,  26, 197,  62,  94, 252, 219, 203, 117,  35,  11,  32,
         57, 177,  33,  88, 237, 149,  56,  87, 174,  20, 125, 136, 171, 168,  68, 175,
         74, 165,  71, 134, 139,  48,  27, 166,  77, 146, 158, 231,  83, 111, 229, 122,
         60, 211, 133, 230, 220, 105,  92,  41,  55,  46, 245,  40, 244, 102, 143,  54,
         65,  25,  63, 161,   1, 216,  80,  73, 209,  76, 132, 187, 208,  89,  18, 169,
        200, 196, 135, 130, 116, 188, 159,  86, 164, 100, 109, 198, 173, 186,   3,  64,
         52, 217, 226, 250, 124, 123,   5, 202,  38, 147, 118, 126, 255,  82,  85, 212,
        207, 206,  59, 227,  47,  16,  58,  17, 182, 189,  28,  42, 223, 183, 170, 213,
        119, 248, 152,   2,  44, 154, 163,  70, 221, 153, 101, 155, 167,  43, 172,   9,
        129,  22,  39, 253,  19,  98, 108, 110,  79, 113, 224, 232, 178, 185, 112, 104,
        218, 246,  97, 228, 251,  34, 242, 193, 238, 210, 144,  12, 191, 179, 162, 241,
         81,  51, 145, 235, 249,  14, 239, 107,  49, 192, 214,  31, 181, 199, 106, 157,
        184,  84, 204, 176, 115, 121,  50,  45, 127,   4, 150, 254, 138, 236, 205,  93,
        222, 114,  67,  29,  24,  72, 243, 141, 128, 195,  78,  66, 215,  61, 156, 180,
    };
}

#endif /* __SLR_distributions__ */
