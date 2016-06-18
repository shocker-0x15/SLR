//
//  voronoi_textures.h
//  SLR
//
//  Created by 渡部 心 on 2016/06/17.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#ifndef voronoi_textures_h
#define voronoi_textures_h

#include "../defines.h"
#include "../references.h"
#include "../Core/textures.h"

namespace SLR {
    class SLR_API VoronoiSpectrumTexture : public SpectrumTexture {
        float m_scale;
    public:
        VoronoiSpectrumTexture(float scale) : m_scale(scale) { }
        
        SampledSpectrum evaluate(const TexCoord2D &tc, const WavelengthSamples &wls) const override;
        RegularConstantContinuous2D* createIBLImportanceMap() const override;
    };
    
    class SLR_API VoronoiNormal3DTexture : public Normal3DTexture {
        float m_scale;
    public:
        VoronoiNormal3DTexture(float scale) : m_scale(scale) { }
        
        Normal3D evaluate(const TexCoord2D &tc) const override;
    };
    
    class SLR_API VoronoiFloatTexture : public FloatTexture {
        float m_scale;
    public:
        VoronoiFloatTexture(float scale) : m_scale(scale) { }
        
        float evaluate(const TexCoord2D &tc) const override;
    };
}

#endif /* voronoi_textures_h */
