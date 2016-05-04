//
//  MicrofacetSurfaceMaterial.h
//  SLR
//
//  Created by 渡部 心 on 2016/05/04.
//  Copyright © 2016年 渡部 心. All rights reserved.
//

#ifndef MicrofacetSurfaceMaterial_h
#define MicrofacetSurfaceMaterial_h

#include "../defines.h"
#include "../references.h"
#include "../Core/surface_material.h"

namespace SLR {
    class SLR_API MicrofacetReflection : public SurfaceMaterial {
        const SVFresnel* m_F;
        const SVMicrofacetDistribution* m_D;
    public:
        MicrofacetReflection(const SVFresnel* F, const SVMicrofacetDistribution* D) :
        m_F(F), m_D(D) {}
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    class SLR_API MicrofacetTransmission : public SurfaceMaterial {
        const SpectrumTexture* m_etaExt;
        const SpectrumTexture* m_etaInt;
        const SVMicrofacetDistribution* m_D;
    public:
        MicrofacetTransmission(const SpectrumTexture* etaExt, const SpectrumTexture* etaInt, const SVMicrofacetDistribution* D) :
        m_etaExt(etaExt), m_etaInt(etaInt), m_D(D) {}
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* MicrofacetSurfaceMaterial_hpp */
