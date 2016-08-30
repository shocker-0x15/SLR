//
//  MicrofacetSurfaceMaterial.cpp
//
//  Created by 渡部 心 on 2016/05/04.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#include "MicrofacetSurfaceMaterial.h"
#include "../BSDFs/MicrofacetBSDF.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/textures.h"

namespace SLR {
    BSDF* MicrofacetReflection::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        SampledSpectrum eta = m_eta->evaluate(surfPt, wls);
        SampledSpectrum k = m_k->evaluate(surfPt, wls);
        MicrofacetDistribution* dist = m_D->getMicrofacetDistribution(surfPt, mem);
        return mem.create<MicrofacetBRDF>(eta, k, dist);
    }
    
    
    
    BSDF* MicrofacetScattering::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        SampledSpectrum etaExt = m_etaExt->evaluate(surfPt, wls);
        SampledSpectrum etaInt = m_etaInt->evaluate(surfPt, wls);
        MicrofacetDistribution* dist = m_D->getMicrofacetDistribution(surfPt, mem);
        return mem.create<MicrofacetBSDF>(etaExt, etaInt, dist);
    }
}
