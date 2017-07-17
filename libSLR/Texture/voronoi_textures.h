//
//  voronoi_textures.h
//
//  Created by 渡部 心 on 2016/06/17.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_voronoi_textures__
#define __SLR_voronoi_textures__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/textures.h"
#include "../Core/distributions.h"

namespace SLR {
    class SLR_API VoronoiSpectrumTexture : public SpectrumTexture {
        const Texture3DMapping* m_mapping;
        WorleyNoise3DGeneratorTemplate<float> m_noiseGen;
        float m_scale;
        float m_brightness;
    public:
        VoronoiSpectrumTexture(const Texture3DMapping* mapping, float scale, float brightness) :
        m_mapping(mapping), m_noiseGen(), m_scale(scale), m_brightness(brightness) { }
        
        SampledSpectrum evaluate(const Point3D &p, const WavelengthSamples &wls) const;
        SampledSpectrum evaluate(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override {
            return evaluate(m_mapping->map(surfPt), wls);
        }
        SampledSpectrum evaluate(const MediumPoint &medPt, const WavelengthSamples &wls) const override {
            return evaluate(m_mapping->map(medPt), wls);
        }
        float evaluateLuminance(const Point3D &p) const;
        float evaluateLuminance(const SurfacePoint &surfPt) const override {
            return evaluateLuminance(m_mapping->map(surfPt) / m_scale);
        }
        float evaluateLuminance(const MediumPoint &medPt) const override {
            return evaluateLuminance(m_mapping->map(medPt) / m_scale);
        }
        const ContinuousDistribution2D* createIBLImportanceMap() const override;
    };
    
    
    
    class SLR_API VoronoiNormalTexture : public NormalTexture {
        const Texture3DMapping* m_mapping;
        WorleyNoise3DGeneratorTemplate<float> m_noiseGen;
        float m_scale;
        float m_cosThetaMax;
    public:
        VoronoiNormalTexture(const Texture3DMapping* mapping, float scale, float thetaMax) :
        m_mapping(mapping), m_noiseGen(), m_scale(scale), m_cosThetaMax(std::cos(thetaMax)) { }
        
        Normal3D evaluate(const Point3D &p) const;
        Normal3D evaluate(const SurfacePoint &surfPt) const override {
            return evaluate(m_mapping->map(surfPt));
        }
        Normal3D evaluate(const MediumPoint &medPt) const override {
            return evaluate(m_mapping->map(medPt));
        }
    };
    
    
    
    class SLR_API VoronoiFloatTexture : public FloatTexture {
        const Texture3DMapping* m_mapping;
        WorleyNoise3DGeneratorTemplate<float> m_noiseGen;
        float m_scale;
        float m_valueScale;
        bool m_flat;
    public:
        VoronoiFloatTexture(const Texture3DMapping* mapping, float scale, float valueScale, bool flat) :
        m_mapping(mapping), m_noiseGen(), m_scale(scale), m_valueScale(valueScale), m_flat(flat) { }
        
        float evaluate(const Point3D &p) const;
        float evaluate(const SurfacePoint &surfPt) const override {
            return evaluate(m_mapping->map(surfPt));
        }
        float evaluate(const MediumPoint &medPt) const override {
            return evaluate(m_mapping->map(medPt));
        }
    };
    
    
    
    class SLR_API WorleyNoiseFloatTexture : public FloatTexture {
        const Texture3DMapping* m_mapping;
        MultiOctaveWorleyNoise3DGeneratorTemplate<float> m_generator;
    public:
        WorleyNoiseFloatTexture(const Texture3DMapping* mapping, uint32_t numOctaves, float initialFrequency, float supValueOrInitialAmplitude, bool supSpecified, float clipValue, 
                                float frequencyMultiplier, float persistence) :  
        m_mapping(mapping), m_generator(numOctaves, initialFrequency, supValueOrInitialAmplitude, supSpecified, clipValue, frequencyMultiplier, persistence) { }
        
        float evaluate(const Point3D &p) const;
        float evaluate(const SurfacePoint &surfPt) const override {
            return evaluate(m_mapping->map(surfPt));
        }
        float evaluate(const MediumPoint &medPt) const override {
            return evaluate(m_mapping->map(medPt));
        }
    };
}

#endif /* __SLR_voronoi_textures__ */
