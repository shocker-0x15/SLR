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
#include "../BasicTypes/Point3D.h"
#include "../BasicTypes/Vector3D.h"

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
    
    
    
    static const uint32_t FNV_OFFSET_BASIS_32 = 2166136261U;
    static const uint64_t FNV_OFFSET_BASIS_64 = 14695981039346656037U;
    
    static const uint32_t FNV_PRIME_32 = 16777619U;
    static const uint64_t FNV_PRIME_64 = 1099511628211LLU;
    
    static inline uint32_t getFNV1Hash32(uint8_t *bytes, size_t length) {
        uint32_t hash = FNV_OFFSET_BASIS_32;
        for (int i = 0; i < length; ++i)
            hash = (FNV_PRIME_32 * hash) ^ (bytes[i]);
        
        return hash;
    }
    
    static inline uint64_t getFNV1Hash64(uint8_t *bytes, size_t length) {
        uint64_t hash = FNV_OFFSET_BASIS_64;
        for (int i = 0; i < length; ++i)
            hash = (FNV_PRIME_64 * hash) ^ (bytes[i]);
        
        return hash;
    }
    
    
    
    // References
    // Improving Noise
    // This code is based on the web site: adrian's soapbox
    // http://flafla2.github.io/2014/08/09/perlinnoise.html
    template <typename RealType>
    class ImprovedPerlinNoise3DGeneratorTemplate {
        static const uint8_t PermutationTable[];
        
        int32_t m_repeat;
        
        static uint8_t hash(int32_t x, int32_t y, int32_t z);
        static RealType gradient(uint32_t hash, RealType xu, RealType yu, RealType zu);
        
    public:
        ImprovedPerlinNoise3DGeneratorTemplate(int32_t repeat) : m_repeat(repeat) {}
        
        RealType evaluate(const Point3DTemplate<RealType> &p) const;
    };
    
    // Reference:
    // Long-Period Hash Functions for Procedural Texturing
    // Permutation table of the hash function of period 739,024 = lcm(11, 13, 16, 17, 19)
    template <typename RealType>
    const uint8_t ImprovedPerlinNoise3DGeneratorTemplate<RealType>::PermutationTable[] = {
        // table 0: 11 numbers
        0, 10, 2, 7, 3, 5, 6, 4, 8, 1, 9, 
        // table 1: 13 numbers
        5, 11, 6, 8, 1, 10, 12, 9, 3, 7, 0, 4, 2, 
        // table 2: 16 numbers = the range of the hash function required by Perlin noise.
        13, 10, 11, 5, 6, 9, 4, 3, 8, 7, 14, 2, 0, 1, 15, 12, 
        // table 3: 17 numbers
        1, 13, 5, 14, 12, 3, 6, 16, 0, 8, 9, 2, 11, 4, 15, 7, 10, 
        // table 4: 19 numbers
        10, 6, 5, 8, 15, 0, 17, 7, 14, 18, 13, 16, 2, 9, 12, 1, 11, 4, 3,  
//        // table 6: 23 numbers
//        20, 21, 4, 5, 0, 18, 14, 2, 6, 22, 10, 17, 3, 7, 8, 16, 19, 11, 9, 13, 1, 15, 12
    };
    
    
    
    template <typename RealType>
    class MultiOctavePerlinNoise3DGeneratorTemplate {
        ImprovedPerlinNoise3DGeneratorTemplate<RealType> m_primaryNoiseGen;
        uint32_t m_numOctaves;
        RealType m_initialFrequency;
        RealType m_initialAmplitude;
        RealType m_frequencyMultiplier;
        RealType m_persistence;
        RealType m_supValue;
        
    public:
        MultiOctavePerlinNoise3DGeneratorTemplate(uint32_t numOctaves, RealType initialFrequency, RealType supValueOrInitialAmplitude, bool supSpecified,  
                                                  RealType frequencyMultiplier, RealType persistence, int32_t repeat) :
        m_primaryNoiseGen(repeat), 
        m_numOctaves(numOctaves), 
        m_initialFrequency(initialFrequency), 
        m_frequencyMultiplier(frequencyMultiplier), m_persistence(persistence) {
            if (supSpecified) {
                RealType amplitude = 1.0f;
                RealType tempSupValue = 0;
                for (int i = 0; i < m_numOctaves; ++i) {                
                    tempSupValue += amplitude;
                    amplitude *= m_persistence;
                }
                m_initialAmplitude = supValueOrInitialAmplitude / tempSupValue;
                m_supValue = supValueOrInitialAmplitude;
            }
            else {
                m_initialAmplitude = supValueOrInitialAmplitude;
                RealType amplitude = m_initialAmplitude;
                m_supValue = 0;
                for (int i = 0; i < m_numOctaves; ++i) {                
                    m_supValue += amplitude;
                    amplitude *= m_persistence;
                }
            }
        }
        
        RealType getSupValue() const {
            return m_supValue;
        }
        
        RealType evaluate(const Point3DTemplate<RealType> &p) const;
    };
    
    
    
    // References:
    // Curl-Noise for Procedural Fluid Flow
    template <typename RealType>
    class CurlNoise3DGeneratorTemplate {
        static constexpr RealType InitialAmplitude = 1.0f;
        static constexpr RealType FrequencyMultiplier = 2.0f;
        static constexpr RealType Persistence = 0.5f;
        
        MultiOctavePerlinNoise3DGeneratorTemplate<RealType> m_noiseGen;
        RealType m_maxFrequency;
        
    public:
        CurlNoise3DGeneratorTemplate(uint32_t numOctaves, RealType initialFrequency) : 
        m_noiseGen(numOctaves, initialFrequency, InitialAmplitude, false, FrequencyMultiplier, Persistence, -1), 
        m_maxFrequency(initialFrequency * std::pow(FrequencyMultiplier, numOctaves - 1)) {
        }
        
        Vector3DTemplate<RealType> evaluate(const Point3DTemplate<RealType> &p) const;
    };
    
    
    
    // References
    // A Cellular Texture Basis Function
    template <typename RealType>
    class WorleyNoise3DGeneratorTemplate {
    public:
        WorleyNoise3DGeneratorTemplate()  { }
        
        void evaluate(const Point3DTemplate<RealType> &p, RealType* closestSqDistance, uint32_t* hashOfClosest, uint32_t* closestFPIdx) const;
    };
    
    
    
    template <typename RealType>
    class MultiOctaveWorleyNoise3DGeneratorTemplate {
        WorleyNoise3DGeneratorTemplate<RealType> m_primaryNoiseGen;
        uint32_t m_numOctaves;
        RealType m_initialFrequency;
        RealType m_initialAmplitude;
        RealType m_frequencyMultiplier;
        RealType m_persistence;
        RealType m_supValue;
        RealType m_clipValue;
        
    public:
        MultiOctaveWorleyNoise3DGeneratorTemplate(uint32_t numOctaves, RealType initialFrequency, RealType supValueOrInitialAmplitude, bool supSpecified, RealType clipValue,  
                                                  RealType frequencyMultiplier, RealType persistence) :
        m_primaryNoiseGen(), 
        m_numOctaves(numOctaves), 
        m_initialFrequency(initialFrequency), 
        m_clipValue(clipValue), 
        m_frequencyMultiplier(frequencyMultiplier), m_persistence(persistence) {
            if (supSpecified) {
                RealType amplitude = 1.0f;
                RealType tempSupValue = 0;
                for (int i = 0; i < m_numOctaves; ++i) {                
                    tempSupValue += amplitude;
                    amplitude *= m_persistence;
                }
                m_initialAmplitude = supValueOrInitialAmplitude / tempSupValue;
                m_supValue = supValueOrInitialAmplitude;
            }
            else {
                m_initialAmplitude = supValueOrInitialAmplitude;
                RealType amplitude = m_initialAmplitude;
                m_supValue = 0;
                for (int i = 0; i < m_numOctaves; ++i) {                
                    m_supValue += amplitude;
                    amplitude *= m_persistence;
                }
            }
        }
        
        RealType getSupValue() const {
            return m_supValue;
        }
        
        RealType evaluate(const Point3DTemplate<RealType> &p) const;
    };
}

#endif /* __SLR_distributions__ */
