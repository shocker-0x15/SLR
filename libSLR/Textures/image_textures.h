//
//  image_textures.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__image_textures__
#define __SLR__image_textures__

#include "../defines.h"
#include "../references.h"
#include "../Core/textures.h"

namespace SLR {
    class SLR_API ImageSpectrumTexture : public SpectrumTexture {
        const TiledImage2D* m_data;
        const Texture2DMapping* m_mapping;
    public:
        ImageSpectrumTexture(const TiledImage2D* image, const Texture2DMapping* mapping) :
        m_data(image), m_mapping(mapping) { }
        
        SampledSpectrum evaluate(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override;
        RegularConstantContinuous2D* createIBLImportanceMap() const override;
    };
    
    class SLR_API ImageNormal3DTexture : public Normal3DTexture {
        const TiledImage2D* m_data;
        const Texture2DMapping* m_mapping;
    public:
        ImageNormal3DTexture(const TiledImage2D* image, const Texture2DMapping* mapping) :
        m_data(image), m_mapping(mapping) { }
        
        Normal3D evaluate(const SurfacePoint &surfPt) const override;
    };
    
    class SLR_API ImageFloatTexture : public FloatTexture {
        const TiledImage2D* m_data;
        const Texture2DMapping* m_mapping;
    public:
        ImageFloatTexture(const TiledImage2D* image, const Texture2DMapping* mapping) :
        m_data(image), m_mapping(mapping) { }
        
        float evaluate(const SurfacePoint &surfPt) const override;
    };
}

#endif /* defined(__SLR__image_textures__) */
