//
//  basic_SurfaceMaterials.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "basic_SurfaceMaterials.h"
#include "../BSDFs/basic_BSDFs.h"
#include "../Memory/ArenaAllocator.h"
#include "textures.h"

BSDF* DiffuseReflection::getBSDF(const SurfacePoint &surfPt, ArenaAllocator &mem, float scale) const {
    BSDF* bsdf = nullptr;
    if (m_sigma) {
        float sigma = m_sigma->evaluate(surfPt.texCoord);
        bsdf = nullptr;
    }
    else {
        bsdf = mem.create<LambertianBRDF>(scale * m_reflectance->evaluate(surfPt.texCoord));
    }
    return bsdf;
}

BSDF* SpecularReflection::getBSDF(const SurfacePoint &surfPt, ArenaAllocator &mem, float scale) const {
    Spectrum coeffR = m_coeffR->evaluate(surfPt.texCoord);
    const Fresnel* fresnel = m_fresnel->getFresnel(surfPt, mem);
    return mem.create<SpecularBRDF>(scale * coeffR, fresnel);
}

BSDF* SpecularTransmission::getBSDF(const SurfacePoint &surfPt, ArenaAllocator &mem, float scale) const {
    Spectrum coeffT = m_coeffT->evaluate(surfPt.texCoord);
    float etaExt = m_etaExt->evaluate(surfPt.texCoord);
    float etaInt = m_etaInt->evaluate(surfPt.texCoord);
    return mem.create<SpecularBTDF>(scale * coeffT, etaExt, etaInt);
}
