//
//  AshikhminShirleyReflection.h
//
//  Created by 渡部 心 on 2015/09/06.
//  Copyright (c) 2015年 渡部 心. All rights reserved.
//

#ifndef __SLR_AshikhminShirleyReflection__
#define __SLR_AshikhminShirleyReflection__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/surface_material.h"

namespace SLR {
    class SLR_API AshikhminShirleyReflectionSurfaceMaterial : public SurfaceMaterial {
        const SpectrumTexture* m_Rs;
        const SpectrumTexture* m_Rd;
        const FloatTexture* m_nu;
        const FloatTexture* m_nv;
    public:
        AshikhminShirleyReflectionSurfaceMaterial(const SpectrumTexture* Rs, const SpectrumTexture* Rd, const FloatTexture* nu, const FloatTexture* nv) :
        m_Rs(Rs), m_Rd(Rd), m_nu(nu), m_nv(nv) {}
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* __SLR_AshikhminReflection__ */
