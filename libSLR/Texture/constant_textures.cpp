//
//  constant_textures.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "constant_textures.h"

namespace SLR {
    const ContinuousDistribution2D* ConstantSpectrumTexture::createIBLImportanceMap() const {
        SLRAssert_NotImplemented();
        return nullptr;
    }
    
    void ConstantSpectrumTexture::generateLuminanceChannel() {
        float XYZ[3];
        m_value->convertToXYZ(XYZ);
        m_luminance = XYZ[1];
    }
}
