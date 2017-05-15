//
//  disney_surface_materials.cpp
//
//  Created by 渡部 心 on 2017/05/14.
//  Copyright (c) 2017年 渡部 心. All rights reserved.
//

#include "disney_surface_materials.h"

#include "../MemoryAllocators/ArenaAllocator.h"
#include "../Core/textures.h"
#include "../BSDF/disney_bsdfs.h"

namespace SLR {
    BSDF* DisneyReflectionSurfaceMaterial::getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale) const {
        SLRAssert_NotImplemented();
        return nullptr;
//        return mem.create<DisneyBRDF>(scale * m_baseColor->evaluate(surfPt, wls), scale * m_baseColor->evaluateLuminance(surfPt), 
//                                      m_subsurface->evaluate(surfPt), m_metallic->evaluate(surfPt), m_specular->evaluate(surfPt), m_specularTint->evaluate(surfPt), 
//                                      m_roughness->evaluate(surfPt), m_anisotropic->evaluate(surfPt), m_sheen->evaluate(surfPt), m_sheenTint->evaluate(surfPt), 
//                                      m_clearCoat->evaluate(surfPt), m_clearCoatGloss->evaluate(surfPt));
    }
}
