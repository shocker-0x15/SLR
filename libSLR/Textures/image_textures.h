//
//  image_textures.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_image_textures__
#define __SLR_image_textures__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/textures.h"

namespace SLR {
    class SLR_API ImageSpectrumTexture : public SpectrumTexture {
        const TiledImage2D* m_data;
        const Texture2DMapping* m_mapping;
    public:
        ImageSpectrumTexture(const TiledImage2D* image, const Texture2DMapping* mapping) :
        m_data(image), m_mapping(mapping) { }
        
        SampledSpectrum evaluate(const Point3D &p, const WavelengthSamples &wls) const;
        SampledSpectrum evaluate(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override {
            return evaluate(m_mapping->map(surfPt), wls);
        }
        SampledSpectrum evaluate(const MediumPoint &medPt, const WavelengthSamples &wls) const override {
            return evaluate(m_mapping->map(medPt), wls);
        }
        RegularConstantContinuousDistribution2D* createIBLImportanceMap() const override;
    };
    
    class SLR_API ImageNormalTexture : public NormalTexture {
        const TiledImage2D* m_data;
        const Texture2DMapping* m_mapping;
    public:
        ImageNormalTexture(const TiledImage2D* image, const Texture2DMapping* mapping) :
        m_data(image), m_mapping(mapping) { }
        
        Normal3D evaluate(const Point3D &p) const ;
        Normal3D evaluate(const SurfacePoint &surfPt) const override {
            return evaluate(m_mapping->map(surfPt));
        }
        Normal3D evaluate(const MediumPoint &medPt) const override {
            return evaluate(m_mapping->map(medPt));
        }
    };
    
    class SLR_API ImageFloatTexture : public FloatTexture {
        const TiledImage2D* m_data;
        const Texture2DMapping* m_mapping;
    public:
        ImageFloatTexture(const TiledImage2D* image, const Texture2DMapping* mapping) :
        m_data(image), m_mapping(mapping) { }
        
        float evaluate(const Point3D &p) const;
        float evaluate(const SurfacePoint &surfPt) const override {
            return evaluate(m_mapping->map(surfPt));
        }
        float evaluate(const MediumPoint &medPt) const override {
            return evaluate(m_mapping->map(medPt));
        }
    };
}

#endif /* __SLR_image_textures__ */
