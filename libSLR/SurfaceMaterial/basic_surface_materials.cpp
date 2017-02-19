//
//  basic_surface_materials.cpp
//
//  Created by 渡部 心 on 2015/09/05.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "basic_surface_materials.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/textures.h"
#include "../BSDF/basic_bsdfs.h"
#include "../BSDF/OrenNayerBRDF.h"

namespace SLR {
    BSDF* DiffuseReflectionSurfaceMaterial::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        BSDF* bsdf = nullptr;
        if (m_sigma) {
            float sigma = m_sigma->evaluate(surfPt);
            bsdf = mem.create<OrenNayerBRDF>(scale * m_reflectance->evaluate(surfPt, wls), sigma);
        }
        else {
            bsdf = mem.create<LambertianBRDF>(scale * m_reflectance->evaluate(surfPt, wls));
        }
        return bsdf;
    }
    
    
    
    BSDF* SpecularReflectionSurfaceMaterial::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        SampledSpectrum coeffR = m_coeffR->evaluate(surfPt, wls);
        SampledSpectrum eta = m_eta->evaluate(surfPt, wls);
        SampledSpectrum k = m_k->evaluate(surfPt, wls);
        return mem.create<SpecularBRDF>(scale * coeffR, eta, k);
    }
    
    
    
    BSDF* SpecularScatteringSurfaceMaterial::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        SampledSpectrum coeff = m_coeff->evaluate(surfPt, wls);
        SampledSpectrum etaExt = m_etaExt->evaluate(surfPt, wls);
        SampledSpectrum etaInt = m_etaInt->evaluate(surfPt, wls);
        return mem.create<SpecularBSDF>(scale * coeff, etaExt, etaInt, !wls.lambdaSelected());
    }
    
    
    
    BSDF* FlippedSurfaceMaterial::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        BSDF* baseBSDF = m_baseMat->getBSDF(surfPt, wls, mem, scale);
        return mem.create<FlippedBSDF>(baseBSDF);
    }
}
