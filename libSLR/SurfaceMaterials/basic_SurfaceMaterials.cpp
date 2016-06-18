//
//  basic_SurfaceMaterials.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "basic_SurfaceMaterials.h"
#include "../BSDFs/basic_BSDFs.h"
#include "../BSDFs/OrenNayerBRDF.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/textures.h"

namespace SLR {
    BSDF* DiffuseReflection::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        BSDF* bsdf = nullptr;
        if (m_sigma) {
            float sigma = m_sigma->evaluate(surfPt.texCoord);
            bsdf = mem.create<OrenNayerBRDF>(scale * m_reflectance->evaluate(surfPt.texCoord, wls), sigma);
        }
        else {
            bsdf = mem.create<LambertianBRDF>(scale * m_reflectance->evaluate(surfPt.texCoord, wls));
        }
        return bsdf;
    }
    
    
    
    BSDF* SpecularReflection::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        SampledSpectrum coeffR = m_coeffR->evaluate(surfPt.texCoord, wls);
        SampledSpectrum eta = m_eta->evaluate(surfPt.texCoord, wls);
        SampledSpectrum k = m_k->evaluate(surfPt.texCoord, wls);
        return mem.create<SpecularBRDF>(scale * coeffR, eta, k);
    }
    
    
    
    BSDF* SpecularScattering::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        SampledSpectrum coeff = m_coeff->evaluate(surfPt.texCoord, wls);
        SampledSpectrum etaExt = m_etaExt->evaluate(surfPt.texCoord, wls);
        SampledSpectrum etaInt = m_etaInt->evaluate(surfPt.texCoord, wls);
        return mem.create<SpecularBSDF>(scale * coeff, etaExt, etaInt, !wls.lambdaSelected());
    }
    
    
    
    BSDF* InverseSurfaceMaterial::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        BSDF* baseBSDF = m_baseMat->getBSDF(surfPt, wls, mem, scale);
        return mem.create<InverseBSDF>(baseBSDF);
    }
}
