//
//  basic_SurfaceMaterials.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "basic_SurfaceMaterials.h"
#include "../BSDFs/basic_BSDFs.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/textures.h"

namespace SLR {
    BSDF* DiffuseReflection::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        BSDF* bsdf = nullptr;
        if (m_sigma) {
            float sigma = m_sigma->evaluate(surfPt.texCoord);
            bsdf = nullptr;
            SLRAssert(false, "Oren-Nayer BSDF is not implemented.");
//            bsdf = mem.create<OrenNayerBRDF>(scale * m_reflectance->evaluate(surfPt.texCoord, wls), sigma);
        }
        else {
            bsdf = mem.create<LambertianBRDF>(scale * m_reflectance->evaluate(surfPt.texCoord, wls));
        }
        return bsdf;
    }
    
    
    
    BSDF* SpecularReflection::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        SampledSpectrum coeffR = m_coeffR->evaluate(surfPt.texCoord, wls);
        const Fresnel* fresnel = m_fresnel->getFresnel(surfPt, wls, mem);
        return mem.create<SpecularBRDF>(scale * coeffR, fresnel);
    }
    
    
    
    BSDF* SpecularTransmission::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        SampledSpectrum coeffT = m_coeffT->evaluate(surfPt.texCoord, wls);
        SampledSpectrum etaExt = m_etaExt->evaluate(surfPt.texCoord, wls);
        SampledSpectrum etaInt = m_etaInt->evaluate(surfPt.texCoord, wls);
        return mem.create<SpecularBTDF>(scale * coeffT, etaExt, etaInt, !wls.lambdaSelected());
    }
    
    
    
    BSDF* InverseSurfaceMaterial::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        BSDF* baseBSDF = m_baseMat->getBSDF(surfPt, wls, mem, scale);
        return mem.create<InverseBSDF>(baseBSDF);
    }
}
