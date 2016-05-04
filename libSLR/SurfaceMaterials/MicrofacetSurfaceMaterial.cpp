//
//  MicrofacetSurfaceMaterial.cpp
//  SLR
//
//  Created by 渡部 心 on 2016/05/04.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#include "MicrofacetSurfaceMaterial.h"
#include "../BSDFs/MicrofacetBSDF.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/textures.h"

namespace SLR {
    BSDF* MicrofacetReflection::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        Fresnel* fresnel = m_F->getFresnel(surfPt, wls, mem);
        MicrofacetDistribution* dist = m_D->getMicrofacetDistribution(surfPt, mem);
        return mem.create<MicrofacetBRDF>(fresnel, dist);
    }
    
    
    
    BSDF* MicrofacetTransmission::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        SampledSpectrum etaExt = m_etaExt->evaluate(surfPt.texCoord, wls);
        SampledSpectrum etaInt = m_etaInt->evaluate(surfPt.texCoord, wls);
        MicrofacetDistribution* dist = m_D->getMicrofacetDistribution(surfPt, mem);
        return mem.create<MicrofacetBTDF>(etaExt, etaInt, !wls.lambdaSelected(), dist);
    }
}
