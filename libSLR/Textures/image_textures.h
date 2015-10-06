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
    class ImageSpectrumTexture : public SpectrumTexture {
        const TiledImage2D* m_data;
    public:
        ImageSpectrumTexture(const TiledImage2D* image) :
        m_data(image) { };
        
        SampledSpectrum evaluate(const TexCoord2D &tc, const WavelengthSamples &wls) const override;
        RegularConstantContinuous2D* createIBLImportanceMap() const override;
    };
    
    class ImageNormal3DTexture : public Normal3DTexture {
        const TiledImage2D* m_data;
    public:
        ImageNormal3DTexture(const TiledImage2D* image) :
        m_data(image) { };
        
        Normal3D evaluate(const TexCoord2D &tc) const override;
    };
    
    class ImageFloatTexture : public FloatTexture {
        const TiledImage2D* m_data;
    public:
        ImageFloatTexture(const TiledImage2D* image) :
        m_data(image) { };
        
        float evaluate(const TexCoord2D &tc) const override;
    };
}

#endif /* defined(__SLR__image_textures__) */
