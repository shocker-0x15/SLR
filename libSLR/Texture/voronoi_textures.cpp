//
//  voronoi_textures.cpp
//
//  Created by 渡部 心 on 2016/06/17.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "voronoi_textures.h"

#include "../Core/distributions.h"
#include "../Core/geometry.h"
#include "../RNG/LinearCongruentialRNG.h"

namespace SLR {
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
    
    static void evaluateVoronoi(const Point3D &p, float* closestDistance, uint32_t* hashOfClosest, uint32_t* closestFPIdx) {
        int32_t iEvalCoord[3];
        iEvalCoord[0] = std::floor(p.x);
        iEvalCoord[1] = std::floor(p.y);
        iEvalCoord[2] = std::floor(p.z);
        
        int32_t rangeBaseX = -1 + std::round(p.x - iEvalCoord[0]);
        int32_t rangeBaseY = -1 + std::round(p.y - iEvalCoord[1]);
        int32_t rangeBaseZ = -1 + std::round(p.z - iEvalCoord[2]);
        
        *closestDistance = INFINITY;
        for (int iz = rangeBaseZ; iz < rangeBaseZ + 2; ++iz) {
            for (int iy = rangeBaseY; iy < rangeBaseY + 2; ++iy) {
                for (int ix = rangeBaseX; ix < rangeBaseX + 2; ++ix) {
                    int32_t iCoord[3] = {iEvalCoord[0] + ix, iEvalCoord[1] + iy, iEvalCoord[2] + iz};
                    uint32_t hash = getFNV1Hash32((uint8_t*)iCoord, sizeof(iCoord));
                    LinearCongruentialRNG rng(hash);
                    
                    uint32_t numFeaturePoints = 1 + std::min(int32_t(8 * rng.getFloat0cTo1o()), 8);
                    for (int i = 0; i < numFeaturePoints; ++i) {
                        Point3D fp = Point3D(iCoord[0] + rng.getFloat0cTo1o(), iCoord[1] + rng.getFloat0cTo1o(), iCoord[2] + rng.getFloat0cTo1o());
                        float dist = distance(p, fp);
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
    
    
    
    SampledSpectrum VoronoiSpectrumTexture::evaluate(const Point3D &p, const WavelengthSamples &wls) const {
        float closestDistance;
        uint32_t hash;
        uint32_t fpIdx;
        evaluateVoronoi(p, &closestDistance, &hash, &fpIdx);
        
        LinearCongruentialRNG rng(hash + fpIdx);
        float rgb[3] = {
            rng.getFloat0cTo1o() * m_brightness,
            rng.getFloat0cTo1o() * m_brightness,
            rng.getFloat0cTo1o() * m_brightness
        };
#ifdef Use_Spectral_Representation
        UpsampledContinuousSpectrum spectrum(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, rgb[0], rgb[1], rgb[2]);
        return spectrum.evaluate(wls);
#else
        return SampledSpectrum(rgb[0], rgb[1], rgb[2]);
#endif
    }
    
    const ContinuousDistribution2D* VoronoiSpectrumTexture::createIBLImportanceMap() const {
        SLRAssert_NotImplemented();
        return nullptr;
    }
    
    Normal3D VoronoiNormalTexture::evaluate(const Point3D &p) const {
        float closestDistance;
        uint32_t hash;
        uint32_t fpIdx;
        evaluateVoronoi(p, &closestDistance, &hash, &fpIdx);
        
        LinearCongruentialRNG rng(hash + fpIdx);
        return uniformSampleCone(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), m_cosThetaMax);
    }
    
    float VoronoiFloatTexture::evaluate(const Point3D &p) const {
        float closestDistance;
        uint32_t hash;
        uint32_t fpIdx;
        evaluateVoronoi(p, &closestDistance, &hash, &fpIdx);
        
        if (m_flat) {
            LinearCongruentialRNG rng(hash + fpIdx);
            return m_valueScale * rng.getFloat0cTo1o();
        }
        else {
            return (closestDistance / (1.414213562 * m_scale)) * m_valueScale;
        }
    }
}
