//
//  constant_textures.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "constant_textures.h"
#include "../BasicTypes/Spectrum.h"

Spectrum ConstantSpectrumTexture::evaluate(const TexCoord2D &tc, const WavelengthSamples &wls) const {
    return m_value->evaluate(wls);
};

RegularConstantContinuous2D* ConstantSpectrumTexture::createIBLImportanceMap() const {
    return nullptr;
}
