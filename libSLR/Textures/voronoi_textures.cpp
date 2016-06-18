//
//  voronoi_textures.cpp
//  SLR
//
//  Created by 渡部 心 on 2016/06/17.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#include "voronoi_textures.h"
#include "geometry.h"
#include "../Core/LinearCongruentialRNG.h"

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
    
    static uint64_t getFNV1Hash64(uint8_t *bytes, size_t length) {
        uint64_t hash = FNV_OFFSET_BASIS_64;
        for (int i = 0; i < length; ++i)
            hash = (FNV_PRIME_64 * hash) ^ (bytes[i]);
        
        return hash;
    }
    
    SampledSpectrum VoronoiSpectrumTexture::evaluate(const TexCoord2D &tc, const WavelengthSamples &wls) const {
        Point3D evalp = Point3D(tc.u, tc.v, 0.0f) / m_scale;
        int32_t iEvalCoord[3];
        iEvalCoord[0] = (int32_t)evalp.x;
        iEvalCoord[1] = (int32_t)evalp.y;
        iEvalCoord[2] = (int32_t)evalp.z;
        
        int32_t rangeBaseX = -1 + std::round(evalp.x - iEvalCoord[0]);
        int32_t rangeBaseY = -1 + std::round(evalp.y - iEvalCoord[1]);
        int32_t rangeBaseZ = -1 + std::round(evalp.z - iEvalCoord[2]);
        
        float closestDistance = INFINITY;
        for (int iz = rangeBaseZ; iz < rangeBaseZ + 2; ++iz) {
            for (int iy = rangeBaseY; iy < rangeBaseY + 2; ++iy) {
                for (int ix = rangeBaseX; ix < rangeBaseX + 2; ++ix) {
                    int32_t iCoord[3] = {iEvalCoord[0] + ix, iEvalCoord[1] + iy, iEvalCoord[2] + iz};
                    LinearCongruentialRNG rng(getFNV1Hash32((uint8_t*)iCoord, sizeof(iCoord)));
                    
                    uint32_t numFeaturePoints = 1 + std::min(int32_t(8 * rng.getFloat0cTo1o()), 8);
                    for (int i = 0; i < numFeaturePoints; ++i) {
                        Point3D fp = Point3D(iCoord[0] + rng.getFloat0cTo1o(), iCoord[1] + rng.getFloat0cTo1o(), iCoord[2] + rng.getFloat0cTo1o());
                        float dist = distance(evalp, fp);
                        if (dist < closestDistance)
                            closestDistance = dist;
                    }
                }
            }
        }
        
        return SampledSpectrum(closestDistance);
    }
    
    RegularConstantContinuous2D* VoronoiSpectrumTexture::createIBLImportanceMap() const {
        SLRAssert_NotImplemented();
        return nullptr;
    }
    
    float VoronoiFloatTexture::evaluate(const TexCoord2D &tc) const {
        return 0.0f;
    }
    
    Normal3D VoronoiNormal3DTexture::evaluate(const TexCoord2D &tc) const {
        return Normal3D();
    }
}
