//
//  constant_textures.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__constant_textures__
#define __SLR__constant_textures__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/textures.h"

namespace SLR {
    class SLR_API ConstantSpectrumTexture : public SpectrumTexture {
        const AssetSpectrum* m_value;
    public:
        ConstantSpectrumTexture(const AssetSpectrum* value) : m_value(value) { }
        
        SampledSpectrum evaluate(const SurfacePoint &surfPt, const WavelengthSamples &wls) const override { return m_value->evaluate(wls); }
        SampledSpectrum evaluate(const MediumPoint &medPt, const WavelengthSamples &wls) const override { return m_value->evaluate(wls); }
        RegularConstantContinuousDistribution2D* createIBLImportanceMap() const override;
    };
    
    
    
    class SLR_API ConstantFloatTexture : public FloatTexture {
        float m_value;
    public:
        ConstantFloatTexture(float value) : m_value(value) { }
        
        float evaluate(const SurfacePoint &surfPt) const override { return m_value; }
        float evaluate(const MediumPoint &medPt) const override { return m_value; }
    };    
}

#endif /* defined(__SLR__constant_textures__) */
