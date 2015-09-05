//
//  AshikhminReflection.cpp
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "AshikhminReflection.h"
#include "../BSDFs/AshikhminBRDF.h"
#include "../Memory/ArenaAllocator.h"
#include "textures.h"

BSDF* AshikhminSpecularReflection::getBSDF(const SurfacePoint &surfPt, ArenaAllocator &mem, float scale) const {
    Spectrum Rs = m_Rs->evaluate(surfPt.texCoord);
    float nu = m_nu->evaluate(surfPt.texCoord);
    float nv = m_nv->evaluate(surfPt.texCoord);
    return mem.create<AshikhminSpecularBRDF>(scale * Rs, nu, nv);
}

BSDF* AshikhminDiffuseReflection::getBSDF(const SurfacePoint &surfPt, ArenaAllocator &mem, float scale) const {
    Spectrum Rs = m_Rs->evaluate(surfPt.texCoord);
    Spectrum Rd = m_Rd->evaluate(surfPt.texCoord);
    return mem.create<AshikhminDiffuseBRDF>(scale * Rs, scale * Rd);
}
