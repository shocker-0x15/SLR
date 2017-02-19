//
//  microfacet_surface_materials.h
//
//  Created by 渡部 心 on 2016/05/04.
//  Copyright (c) 2016年 渡部 心. All rights reserved.
//

#ifndef __SLR_microfacet_surface_materials__
#define __SLR_microfacet_surface_materials__

#include "../defines.h"
#include "../declarations.h"
#include "../Core/surface_material.h"

namespace SLR {
    class SLR_API MicrofacetReflectionSurfaceMaterial : public SurfaceMaterial {
        const SpectrumTexture* m_eta;
        const SpectrumTexture* m_k;
        const SVMicrofacetDistribution* m_D;
    public:
        MicrofacetReflectionSurfaceMaterial(const SpectrumTexture* eta, const SpectrumTexture* k, const SVMicrofacetDistribution* D) :
        m_eta(eta), m_k(k), m_D(D) {}
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
    
    
    
    class SLR_API MicrofacetScatteringSurfaceMaterial : public SurfaceMaterial {
        const SpectrumTexture* m_etaExt;
        const SpectrumTexture* m_etaInt;
        const SVMicrofacetDistribution* m_D;
    public:
        MicrofacetScatteringSurfaceMaterial(const SpectrumTexture* etaExt, const SpectrumTexture* etaInt, const SVMicrofacetDistribution* D) :
        m_etaExt(etaExt), m_etaInt(etaInt), m_D(D) {}
        
        BSDF* getBSDF(const SurfacePoint &surfPt, const WavelengthSamples &wls, ArenaAllocator &mem, float scale = 1.0f) const override;
    };
}

#endif /* __SLR_microfacet_surface_materials__ */
