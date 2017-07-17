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
    SampledSpectrum VoronoiSpectrumTexture::evaluate(const Point3D &p, const WavelengthSamples &wls) const {
        float closestDistance;
        uint32_t hash;
        uint32_t fpIdx;
        m_noiseGen.evaluate(p / m_scale, &closestDistance, &hash, &fpIdx);
        
        LinearCongruentialRNG rng(hash + fpIdx);
        float rgb[3] = {
            rng.getFloat0cTo1o() * m_brightness,
            rng.getFloat0cTo1o() * m_brightness,
            rng.getFloat0cTo1o() * m_brightness
        };
#ifdef SLR_Use_Spectral_Representation
        UpsampledContinuousSpectrum spectrum(SpectrumType::Reflectance, ColorSpace::sRGB_NonLinear, rgb[0], rgb[1], rgb[2]);
        return spectrum.evaluate(wls);
#else
        return SampledSpectrum(rgb[0], rgb[1], rgb[2]);
#endif
    }
    
    float VoronoiSpectrumTexture::evaluateLuminance(const Point3D &p) const {
        float closestDistance;
        uint32_t hash;
        uint32_t fpIdx;
        m_noiseGen.evaluate(p / m_scale, &closestDistance, &hash, &fpIdx);
        
        LinearCongruentialRNG rng(hash + fpIdx);
        float rgb[3] = {
            rng.getFloat0cTo1o() * m_brightness,
            rng.getFloat0cTo1o() * m_brightness,
            rng.getFloat0cTo1o() * m_brightness
        };
        
        return sRGB_to_Luminance(rgb[0], rgb[1], rgb[2]);
    }
    
    const ContinuousDistribution2D* VoronoiSpectrumTexture::createIBLImportanceMap() const {
        SLRAssert_NotImplemented();
        return nullptr;
    }
    
    Normal3D VoronoiNormalTexture::evaluate(const Point3D &p) const {
        float closestDistance;
        uint32_t hash;
        uint32_t fpIdx;
        m_noiseGen.evaluate(p / m_scale, &closestDistance, &hash, &fpIdx);
        
        LinearCongruentialRNG rng(hash + fpIdx);
        return uniformSampleCone(rng.getFloat0cTo1o(), rng.getFloat0cTo1o(), m_cosThetaMax);
    }
    
    float VoronoiFloatTexture::evaluate(const Point3D &p) const {
        float closestDistance;
        uint32_t hash;
        uint32_t fpIdx;
        m_noiseGen.evaluate(p, &closestDistance, &hash, &fpIdx);
        
        if (m_flat) {
            LinearCongruentialRNG rng(hash + fpIdx);
            return m_valueScale * rng.getFloat0cTo1o();
        }
        else {
            return (closestDistance / (1.414213562 * m_scale)) * m_valueScale;
        }
    }
    
    
    
    float WorleyNoiseFloatTexture::evaluate(const Point3D &p) const {
        return m_generator.evaluate(p);
    }
}
