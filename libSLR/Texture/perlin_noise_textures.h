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
    class SLR_API PerlinNoiseFloatTexture : public FloatTexture {
        const Texture3DMapping* m_mapping;
        MultiOctaveImprovedPerlinNoise3DGenerator<float> m_generator;
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
