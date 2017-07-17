//
//  perlin_noise_textures.h
//
//  Created by 渡部 心 on 2017/03/17.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#ifndef __SLR_perlin_noise_textures__
#define __SLR_perlin_noise_textures__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/textures.h"
#include "../Core/distributions.h"

namespace SLR {
    class SLR_API PerlinNoiseNormalTexture : public NormalTexture {
        const Texture3DMapping* m_mapping;
        float m_thetaMax;
        MultiOctavePerlinNoise3DGeneratorTemplate<float> m_generator[2];
    public:
        PerlinNoiseNormalTexture(const Texture3DMapping* mapping, float thetaMax, 
                                 uint32_t numOctaves, float initialFrequencyPhi, float initialFrequencyTheta, 
                                 float frequencyMultiplier, float persistence, int32_t repeat) :  
        m_mapping(mapping), m_thetaMax(thetaMax), 
        m_generator{{numOctaves, initialFrequencyPhi, 1.0f, true, frequencyMultiplier, persistence, repeat}, 
                    {numOctaves, initialFrequencyTheta, 1.0f, true, frequencyMultiplier, persistence, repeat}} { }
        
        Normal3D evaluate(const Point3D &p) const;
        Normal3D evaluate(const SurfacePoint &surfPt) const override {
            return evaluate(m_mapping->map(surfPt));
        }
        Normal3D evaluate(const MediumPoint &medPt) const override {
            return evaluate(m_mapping->map(medPt));
        }
    };
    
    
    
    class SLR_API PerlinNoiseFloatTexture : public FloatTexture {
        const Texture3DMapping* m_mapping;
        MultiOctavePerlinNoise3DGeneratorTemplate<float> m_generator;
    public:
        PerlinNoiseFloatTexture(const Texture3DMapping* mapping, uint32_t numOctaves, float initialFrequency, float supValueOrInitialAmplitude, bool supSpecified, 
                                float frequencyMultiplier, float persistence, int32_t repeat) :  
        m_mapping(mapping), m_generator(numOctaves, initialFrequency, supValueOrInitialAmplitude, supSpecified, frequencyMultiplier, persistence, repeat) { }
        
        float evaluate(const Point3D &p) const;
        float evaluate(const SurfacePoint &surfPt) const override {
            return evaluate(m_mapping->map(surfPt));
        }
        float evaluate(const MediumPoint &medPt) const override {
            return evaluate(m_mapping->map(medPt));
        }
    };
}

#endif /* __SLR_perlin_noise_textures__ */
