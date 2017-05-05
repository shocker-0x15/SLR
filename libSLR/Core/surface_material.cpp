//
//  surface_material.cpp
//
//  Created by 渡部 心 on 2015/04/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "surface_material.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/directional_distribution_functions.h"
#include "../Core/textures.h"

namespace SLR {
    Fresnel* SVFresnelNoOp::getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return mem.create<FresnelNoOp>();
    }
    
    Fresnel* SVFresnelConductor::getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return mem.create<FresnelConductor>(m_eta->evaluate(surfPt, wls), m_k->evaluate(surfPt, wls));
    }
    
    Fresnel* SVFresnelDielectric::getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return mem.create<FresnelDielectric>(m_etaExt->evaluate(surfPt, wls), m_etaInt->evaluate(surfPt, wls));
    }
    
    
    
    MicrofacetDistribution* GGXSVMicrofacetDistribution::getMicrofacetDistribution(const SurfacePoint &surfPt, ArenaAllocator &mem) const {
        return mem.create<GGXMicrofacetDistribution>(m_alpha_gx->evaluate(surfPt.getTextureCoordinate()), 
                                                     m_alpha_gy->evaluate(surfPt.getTextureCoordinate()));
    }
}
