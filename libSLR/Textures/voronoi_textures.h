//
//  voronoi_textures.h
//
//  Created by 渡部 心 on 2016/06/17.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef voronoi_textures_h
#define voronoi_textures_h

#include "../defines.h"
#include "../references.h"
#include "../Core/textures.h"

namespace SLR {
    class SLR_API VoronoiSpectrumTexture : public SpectrumTexture {
        const Texture3DMapping* m_mapping;
        float m_scale;
        float m_brightness;
    public:
        VoronoiSpectrumTexture(const Texture3DMapping* mapping, float scale, float brightness) :
        m_mapping(mapping), m_scale(scale), m_brightness(brightness) { }
        
        SampledSpectrum evaluate(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override;
        RegularConstantContinuous2D* createIBLImportanceMap() const override;
    };
    
    class SLR_API VoronoiNormal3DTexture : public Normal3DTexture {
        const Texture3DMapping* m_mapping;
        float m_scale;
        float m_cosThetaMax;
    public:
        VoronoiNormal3DTexture(const Texture3DMapping* mapping, float scale, float thetaMax) :
        m_mapping(mapping), m_scale(scale), m_cosThetaMax(std::cos(thetaMax)) { }
        
        Normal3D evaluate(const SurfacePoint &surfPt) const override;
    };
    
    class SLR_API VoronoiFloatTexture : public FloatTexture {
        const Texture3DMapping* m_mapping;
        float m_scale;
        float m_valueScale;
        bool m_flat;
    public:
        VoronoiFloatTexture(const Texture3DMapping* mapping, float scale, float valueScale, bool flat) :
        m_mapping(mapping), m_scale(scale), m_valueScale(valueScale), m_flat(flat) { }
        
        float evaluate(const SurfacePoint &surfPt) const override;
    };
}

#endif /* voronoi_textures_h */
