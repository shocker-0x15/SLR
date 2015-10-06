//
//  surface_material.cpp
//
//  Created by 渡部 心 on 2015/04/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#include "surface_material.h"
#include "../Memory/ArenaAllocator.h"
#include "../Core/directional_distribution_functions.h"
#include "../Core/textures.h"

namespace SLR {
    Fresnel* SpatialFresnelNoOp::getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return mem.create<FresnelNoOp>();
    }
    
    Fresnel* SpatialFresnelConductor::getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return mem.create<FresnelConductor>(m_eta->evaluate(surfPt.texCoord, wls), m_k->evaluate(surfPt.texCoord, wls));
    }
    
    Fresnel* SpatialFresnelDielectric::getFresnel(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem) const {
        return mem.create<FresnelDielectric>(m_etaExt->evaluate(surfPt.texCoord, wls), m_etaInt->evaluate(surfPt.texCoord, wls));
    }
    
    BSDF* EmitterSurfaceMaterial::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        return m_mat->getBSDF(surfPt, wls, mem);
    }
    
    SampledSpectrum EmitterSurfaceMaterial::emittance(const SurfacePoint &surfPt, const WavelengthSamples &wls) const {
        return m_emit->emittance(surfPt, wls);
    }
    
    EDF* EmitterSurfaceMaterial::getEDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        return m_emit->getEDF(surfPt, wls, mem);
    }
}
