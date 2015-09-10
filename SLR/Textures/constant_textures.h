//
//  constant_textures.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR__constant_textures__
#define __SLR__constant_textures__

#include "../defines.h"
#include "../references.h"
#include "../Core/textures.h"

class ConstantSpectrumTexture : public SpectrumTexture {
    Spectrum m_value;
public:
    ConstantSpectrumTexture(const Spectrum &value) : m_value(value) { };
    
    Spectrum evaluate(const TexCoord2D &tc, const WavelengthSamples &wls) const override { return m_value; };
    RegularConstantContinuous2D* createIBLImportanceMap() const override;
};

class ConstantFloatTexture : public FloatTexture {
    float m_value;
public:
    ConstantFloatTexture(const float &value) : m_value(value) { };
    
    float evaluate(const TexCoord2D &tc) const override { return m_value; };
};

#endif /* defined(__SLR__constant_textures__) */
