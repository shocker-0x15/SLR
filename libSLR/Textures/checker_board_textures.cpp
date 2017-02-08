//
//  checker_board_textures.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "checker_board_textures.h"

namespace SLR {
    RegularConstantContinuous2D* CheckerBoardSpectrumTexture::createIBLImportanceMap() const {
        SLRAssert_NotImplemented();
        return nullptr;
    }
    
    Normal3D CheckerBoardNormal3DTexture::evaluate(const Point3D &p) const {
        float halfWidth = m_stepWidth * 0.5f;
        
        float uComp = 0.0f;
        float absWrapU = std::fmod(std::fabs(p.x), 1.0f);
        if (absWrapU < halfWidth * 0.5f || absWrapU > 1.0f - halfWidth * 0.5f)
            uComp = 1.0f;
        else if (absWrapU > 0.5f - halfWidth * 0.5f && absWrapU < 0.5f + halfWidth * 0.5f)
            uComp = -1.0f;
        
        float vComp = 0.0f;
        float absWrapV = std::fmod(std::fabs(p.y), 1.0f);
        if (absWrapV < halfWidth * 0.5f || absWrapV > 1.0f - halfWidth * 0.5f)
            vComp = 1.0f;
        else if (absWrapV > 0.5f - halfWidth * 0.5f && absWrapV < 0.5f + halfWidth * 0.5f)
            vComp = -1.0f;
        
        if (absWrapV > 0.5f)
            uComp *= -1;
        if (absWrapU > 0.5f)
            vComp *= -1;
        if (m_reverse) {
            uComp *= -1;
            vComp *= -1;
        }
        
        return normalize(Normal3D(uComp, vComp, 1.0f));
    }    
}
